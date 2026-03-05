"""
Executor for extracting raw signal using Uncalled4 read-level mapping (self mapping).
Combines BAM alignments with Uncalled4 eventalign TSV outputs using a high-performance 
dual-table SQLite database and Producer-Consumer multiprocessing.
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

def worker_extractor(task_queue, result_queue, db_path, bed_target_sites):
    """
    Consumer worker function. Pulls batches of read alignments from the BAM queue, 
    retrieves the precise signal boundaries from the TSV-derived SQLite table,
    extracts the signal from POD5, and writes the output.
    """
    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    except sqlite3.OperationalError:
        conn = sqlite3.connect(db_path)
        
    c = conn.cursor()
    
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
                
                # 1. Get POD5 path
                c.execute("SELECT path FROM pod5_map WHERE read_name = ?", (read_id,))
                res = c.fetchone()
                if not res: 
                    continue
                pod5_path = res[0]
                
                # 2. Recover precise signal map from TSV data (Usando la logica originale)
                c.execute("SELECT query_pos, start_idx, end_idx, kmer FROM signal_map WHERE read_name = ?", (read_id,))
                rows = c.fetchall()
                if not rows: 
                    continue
                    
                # Map TSV position to (start, end, kmer)
                local_signal_map = {r[0]: (r[1], r[2], r[3]) for r in rows}
                
                # 3. Load full read signal
                try:
                    with pod5.Reader(pod5_path) as reader:
                        selection = next(reader.reads(selection=[read_id]))
                        full_signal = signal_math.get_picoampere_signal(selection)
                except Exception:
                    continue
                
                # 4. Extract aligned segments
                signal_len = len(full_signal)
                for q_pos, r_pos in aligned_pairs:
                    # Optional BED filter
                    if bed_target_sites and not bed_target_sites.contains(chrom, r_pos):
                        continue
                    
                    # In Uncalled4 self-mode, the TSV 'position' maps to the read query position (q_pos)
                    signal_info = local_signal_map.get(q_pos)
                    if not signal_info: 
                        continue
                    
                    try:
                        start_idx = int(signal_info[0])
                        end_idx = int(signal_info[1])
                        read_kmer = str(signal_info[2])
                    except (ValueError, TypeError):
                        continue
                    
                    if start_idx < 0 or end_idx > signal_len or start_idx >= end_idx: 
                        continue
                    
                    chunk_pA = full_signal[start_idx:end_idx]
                    sig_str = ",".join([f"{x:.1f}" for x in chunk_pA])
                    
                    # format: read_id, contig, position, ref_kmer(.), read_kmer, samples
                    buffer.append(f"{read_id}\t{chrom}\t{r_pos}\t.\t{read_kmer}\t{sig_str}")
                
            if buffer:
                out_f.write("\n".join(buffer) + "\n")
    
    conn.close()
    result_queue.put(temp_filename)

def setup_databases(args, db_path):
    """
    Creates a dual-table SQLite database indexing both POD5 files and the Uncalled4 TSV.
    """
    logger.info("Setting up SQLite mapping database...")
    if os.path.exists(db_path): 
        os.remove(db_path)
        
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    
    # Performance pragmas for massive inserts
    c.execute("PRAGMA synchronous = OFF")
    c.execute("PRAGMA journal_mode = MEMORY")
    
    # Corretto lo schema per riflettere l'estrazione basata su query_pos
    c.execute("CREATE TABLE signal_map (read_name TEXT, query_pos INTEGER, start_idx INTEGER, end_idx INTEGER, kmer TEXT)")
    c.execute("CREATE TABLE pod5_map (read_name TEXT PRIMARY KEY, path TEXT)")
    
    # 1. Load TSV Signal Map in chunks
    logger.info(f"Loading TSV Signal Map from {args.self_tsv}...")
    try:
        chunk_iter = pd.read_csv(args.self_tsv, sep="\t", chunksize=config.EVENTALIGN_CHUNK_SIZE)
        for chunk in chunk_iter:
            r_col = 'read_name' if 'read_name' in chunk.columns else 'read_id'
            data = list(zip(chunk[r_col], chunk['position'], chunk['start_idx'], chunk['end_idx'], chunk['reference_kmer']))
            c.executemany("INSERT INTO signal_map VALUES (?, ?, ?, ?, ?)", data)
        conn.commit()
    except Exception as e:
        logger.error(f"Failed to load TSV: {e}")
        conn.close()
        raise e

    # 2. Use base module to index POD5 files into the same DB
    conn.close()
    base.index_pod5_files(args.pod5, db_path)
    
    # 3. Create index for fast retrieval
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute("CREATE INDEX idx_signal ON signal_map (read_name)")
    conn.commit()
    conn.close()

def run(args):
    """
    Main entry point for the 'unc-self' extraction subcommand.
    """
    db_path = "uncalled_self_map.db"
    
    # 0. BED loading
    bed_target_sites = bed_handling.load_bed_target_sites(getattr(args, 'bed', None))

    # 1. Setup DB
    setup_databases(args, db_path)

    # 2. Setup Multiprocessing Queues
    manager = multiprocessing.Manager()
    task_queue = manager.Queue(maxsize=config.QUEUE_MAX_SIZE)
    result_queue = manager.Queue()
    
    # 3. Start Consumer Workers
    logger.info(f"Starting {args.threads} worker processes...")
    workers = []
    for _ in range(args.threads):
        p = multiprocessing.Process(
            target=worker_extractor, 
            args=(task_queue, result_queue, db_path, bed_target_sites)
        )
        p.start()
        workers.append(p)

    # 4. Producer: Read BAM
    logger.info("Reading BAM (Producer) and dispatching batches...")
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
    
    # Send poison pills
    for _ in range(args.threads):
        task_queue.put(None)
        
    logger.info("BAM reading completed. Waiting for workers to finish...")

    # 5. Collect Results
    temp_files = [result_queue.get() for _ in range(args.threads)]
    for p in workers:
        p.join()

    # 6. Merge Files
    logger.info("Merging worker outputs...")
    with open(args.output, 'w') as outfile:
        outfile.write("read_id\tcontig\tposition\tref_kmer\tread_kmer\tsamples_pA\n")
        for tf in temp_files:
            if os.path.exists(tf):
                with open(tf, 'r') as infile:
                    shutil.copyfileobj(infile, outfile)
                os.remove(tf)

    if os.path.exists(db_path): 
        os.remove(db_path)
        
    logger.info(f"Uncalled4 self extraction completed successfully: {args.output}")
