"""
Executor for extracting raw signal using F5C event alignment mapping.
Implements a Producer-Consumer architecture using multiprocessing queues,
with concurrent reference FASTA querying for K-mer extraction and precise
signal boundaries from f5c eventalign TSV.
"""
import os
import pysam
import sqlite3
import shutil
import multiprocessing
import time
import logging
import pandas as pd

from core import config
from core import base
from core import signal_math
from core import bed_handling

logger = logging.getLogger(__name__)

def worker_extractor(task_queue, result_queue, db_path, ref_fasta_path, bed_target_sites):
    """
    Consumer worker function. Pulls batches of read alignments, 
    retrieves precise signal boundaries from the SQLite TSV map, 
    extracts the raw signal from POD5, and writes the output.
    """
    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    except sqlite3.OperationalError:
        conn = sqlite3.connect(db_path)
        
    c = conn.cursor()
    ref_fasta = pysam.FastaFile(ref_fasta_path) if ref_fasta_path else None
    
    proc_name = multiprocessing.current_process().name
    temp_filename = f"temp_worker_{proc_name}_{int(time.time())}.tsv"
    
    import pod5
    
    with open(temp_filename, 'w') as out_f:
        while True:
            batch = task_queue.get()
            if batch is None: 
                break
                
            buffer = []
            for read_data in batch:
                read_id = read_data['id']
                chrom = read_data['chrom']
                aligned_pairs = read_data['pairs'] 
                query_sequence = read_data.get('seq', None)
                
                # 1. Get POD5 path
                c.execute("SELECT path FROM pod5_map WHERE read_name = ?", (read_id,))
                res = c.fetchone()
                if not res: 
                    continue
                pod5_path = res[0]
                
                # 2. Get precise boundaries from f5c TSV
                c.execute("SELECT kmer_idx, start, end FROM signal_map WHERE read_name = ?", (read_id,))
                rows = c.fetchall()
                if not rows: 
                    continue
                
                # Map q_pos to (start, end)
                local_signal_map = {r[0]: (r[1], r[2]) for r in rows}
                
                # 3. Load full signal
                try:
                    with pod5.Reader(pod5_path) as reader:
                        selection = next(reader.reads(selection=[read_id]))
                        full_signal = signal_math.get_picoampere_signal(selection)
                except Exception:
                    continue
                
                signal_len = len(full_signal)

                # 4. Extract based on BAM alignments
                for q_pos, r_pos in aligned_pairs:
                    if q_pos is None or r_pos is None: 
                        continue

                    # --- FILTRO BED ---
                    if bed_target_sites and not bed_target_sites.contains(chrom, r_pos):
                        continue
                    # ------------------
                    
                    signal_info = local_signal_map.get(q_pos)
                    if not signal_info: 
                        continue
                    
                    try:
                        start_idx = int(signal_info[0])
                        end_idx = int(signal_info[1])
                    except (ValueError, TypeError):
                        continue
                    
                    if start_idx < 0 or end_idx > signal_len or start_idx >= end_idx: 
                        continue
                    
                    # Estrazione K-mer
                    ref_kmer = "NNNNN"
                    if ref_fasta:
                        ref_kmer = signal_math.get_kmer_from_fasta(ref_fasta, chrom, r_pos)
                    
                    read_kmer = "NNNNN"
                    if query_sequence and 2 <= q_pos < len(query_sequence) - 2:
                        read_kmer = query_sequence[q_pos-2 : q_pos+3]

                    # Slice del segnale
                    chunk_pA = full_signal[start_idx:end_idx]
                    sig_str = ",".join([f"{x:.1f}" for x in chunk_pA])
                    
                    buffer.append(f"{read_id}\t{chrom}\t{r_pos}\t{ref_kmer}\t{read_kmer}\t{sig_str}")
            
            if buffer:
                out_f.write("\n".join(buffer) + "\n")
    
    if ref_fasta: 
        ref_fasta.close()
    conn.close()
    
    result_queue.put(temp_filename)

def setup_databases(args, db_path):
    """
    Creates the dual-table SQLite database indexing POD5 and f5c TSV.
    """
    logger.info("Setting up SQLite mapping database...")
    if os.path.exists(db_path): 
        os.remove(db_path)
        
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute("PRAGMA synchronous = OFF")
    c.execute("PRAGMA journal_mode = MEMORY")
    c.execute("CREATE TABLE signal_map (read_name TEXT, kmer_idx INTEGER, start INTEGER, end INTEGER)")
    c.execute("CREATE TABLE pod5_map (read_name TEXT PRIMARY KEY, path TEXT)")
    
    # 1. Load f5c TSV (dalla tua logica originale robusta)
    logger.info(f"Loading f5c TSV Map from {args.self_tsv}...")
    try:
        cols = ["read_id", "kmer_idx", "start_raw_idx", "end_raw_idx"]
        chunk_iter = pd.read_csv(args.self_tsv, sep="\t", chunksize=config.EVENTALIGN_CHUNK_SIZE, usecols=cols, dtype=str)
        
        for chunk in chunk_iter:
            chunk['kmer_idx'] = pd.to_numeric(chunk['kmer_idx'], errors='coerce')
            chunk['start_raw_idx'] = pd.to_numeric(chunk['start_raw_idx'], errors='coerce')
            chunk['end_raw_idx'] = pd.to_numeric(chunk['end_raw_idx'], errors='coerce')
            chunk.dropna(inplace=True)
            
            data = list(zip(
                chunk['read_id'], 
                chunk['kmer_idx'].astype(int), 
                chunk['start_raw_idx'].astype(int), 
                chunk['end_raw_idx'].astype(int)
            ))
            c.executemany("INSERT INTO signal_map VALUES (?, ?, ?, ?)", data)
        conn.commit()
    except Exception as e:
        logger.error(f"Loading TSV failed: {e}")
        conn.close()
        raise e

    # 2. Index POD5 using base module
    conn.close()
    base.index_pod5_files(args.pod5, db_path)
    
    # 3. Create Index
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute("CREATE INDEX idx_signal ON signal_map (read_name)")
    conn.commit()
    conn.close()

def run(args):
    """
    Main entry point for the 'f5c-self' extraction subcommand.
    """
    db_path = "f5c_self_map.db"
    
    # 0. Load BED target sites
    bed_target_sites = bed_handling.load_bed_target_sites(getattr(args, 'bed', None))
    
    # 1. Setup DB
    setup_databases(args, db_path)
    
    # 2. Setup Multiprocessing
    manager = multiprocessing.Manager()
    task_queue = manager.Queue(maxsize=config.QUEUE_MAX_SIZE)
    result_queue = manager.Queue()
    
    logger.info(f"Starting {args.threads} worker processes for F5C extraction...")
    workers = []
    for _ in range(args.threads):
        p = multiprocessing.Process(
            target=worker_extractor, 
            args=(task_queue, result_queue, db_path, getattr(args, 'fasta', None), bed_target_sites)
        )
        p.start()
        workers.append(p)
        
    # 3. Producer: Read BAM
    logger.info("Reading BAM file and dispatching read batches...")
    bam = pysam.AlignmentFile(args.bam, "rb")
    current_batch = []
    count = 0
    
    for read in bam:
        if read.is_unmapped or read.is_secondary: 
            continue
            
        read_data = {
            'id': read.query_name,
            'chrom': read.reference_name,
            'pairs': read.get_aligned_pairs(matches_only=True),
            'seq': read.query_sequence  # We NEED this for read_kmer!
        }
        current_batch.append(read_data)
        
        if len(current_batch) >= config.READ_BATCH_SIZE:
            task_queue.put(current_batch)
            current_batch = []
            count += config.READ_BATCH_SIZE
            if count % 5000 == 0:
                logger.info(f"Dispatched {count} reads to workers...")
                
    if current_batch: 
        task_queue.put(current_batch)
        
    bam.close()
    
    for _ in range(args.threads):
        task_queue.put(None)
         
    logger.info("BAM parsing completed. Waiting for workers to finish...")
    
    # 4. Collect and Merge
    temp_files = [result_queue.get() for _ in range(args.threads)]
    for p in workers: 
        p.join()
    
    logger.info("Merging temporary worker files into final output...")
    with open(args.output, 'w') as outfile:
        outfile.write("read_id\tcontig\tposition\tref_kmer\tread_kmer\tsamples_pA\n")
        for tf in temp_files:
            if os.path.exists(tf):
                with open(tf, 'r') as infile:
                    shutil.copyfileobj(infile, outfile)
                os.remove(tf)
                  
    if os.path.exists(db_path): 
        os.remove(db_path)
        
    logger.info("F5C self extraction completed successfully.")                except Exception as e:
                    continue
                
                # 3. MOCK EXTRACTION LOGIC FOR F5C
                # Simulating extraction per pair using a fixed slice size of 12 samples.
                sig_idx = 0
                for q_pos, r_pos in aligned_pairs:
                    raw_slice = full_signal[sig_idx : sig_idx + 12]
                    
                    # Extract reference K-mer if the position is mapped
                    ref_kmer = "N" * config.KMER_LENGTH
                    if ref_fasta and r_pos is not None:
                        ref_kmer = signal_math.get_kmer_from_fasta(ref_fasta, chrom, r_pos)
                         
                    if len(raw_slice) > 0:
                        sig_str = ",".join([f"{x:.1f}" for x in raw_slice])
                        # Write: read_id, contig, position, ref_kmer, read_kmer(N), samples
                        out_f.write(f"{read_id}\t{chrom}\t{r_pos}\t{ref_kmer}\tN\t{sig_str}\n")
                        
                    sig_idx += 12
                     
    conn.close()
    if ref_fasta:
        ref_fasta.close()
        
    # Send the temporary filename back to the main process
    result_queue.put(temp_filename)

def run(args):
    """
    Main entry point for the 'f5c-self' extraction subcommand.
    """
    db_path = "f5c_self_map.db"
    
    # 1. Indexing
    base.index_pod5_files(args.pod5, db_path)
    
    # 2. Setup Multiprocessing Queues
    manager = multiprocessing.Manager()
    task_queue = manager.Queue(maxsize=config.QUEUE_MAX_SIZE)
    result_queue = manager.Queue()
    
    # 3. Start Consumer Workers
    logger.info(f"Starting {args.threads} worker processes for F5C extraction...")
    workers = []
    for i in range(args.threads):
        p = multiprocessing.Process(
            target=worker_extractor, 
            args=(task_queue, result_queue, db_path, args.fasta)
        )
        p.start()
        workers.append(p)
        
    # 4. Producer: Read BAM and assign tasks
    logger.info("Reading F5C BAM file and dispatching read batches...")
    bam = pysam.AlignmentFile(args.bam, "rb")
    current_batch = []
    count = 0
    
    for read in bam:
        if read.is_unmapped or read.is_secondary: 
            continue
            
        read_data = {
            'id': read.query_name,
            'chrom': read.reference_name,
            'pairs': read.get_aligned_pairs(matches_only=True)
            # The original code extracted 'seq' here, but it wasn't used downstream.
            # Keeping the payload lightweight.
        }
        current_batch.append(read_data)
        
        if len(current_batch) >= config.READ_BATCH_SIZE:
            task_queue.put(current_batch)
            current_batch = []
            count += config.READ_BATCH_SIZE
            if count % 5000 == 0:
                logger.info(f"Dispatched {count} reads to workers...")
                
    # Flush remaining reads
    if current_batch: 
        task_queue.put(current_batch)
        
    bam.close()
    
    # Send poison pills
    for _ in range(args.threads):
        task_queue.put(None)
         
    logger.info("BAM parsing completed. Waiting for workers to finish extraction...")
    
    # 5. Collect Results
    temp_files = [result_queue.get() for _ in range(args.threads)]
    
    for p in workers: 
        p.join()
    
    # 6. Merge Temporary Files
    logger.info("Merging temporary worker files into final output...")
    with open(args.output, 'w') as outfile:
        outfile.write("read_id\tcontig\tposition\tref_kmer\tread_kmer\tsamples_pA\n")
        
        for tf in temp_files:
            if os.path.exists(tf):
                with open(tf, 'r') as infile:
                    shutil.copyfileobj(infile, outfile)
                os.remove(tf)
                  
    # Cleanup SQLite DB
    if os.path.exists(db_path): 
        os.remove(db_path)
        
    logger.info("F5C self extraction completed successfully.")
