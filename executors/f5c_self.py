"""
Executor for extracting raw signal using F5C event alignment mapping.
Implements a Producer-Consumer architecture using multiprocessing queues,
with the addition of concurrent reference FASTA querying for K-mer extraction.
"""
import os
import pysam
import sqlite3
import shutil
import multiprocessing
import time
import logging

from core import config
from core import base
from core import signal_math

logger = logging.getLogger(__name__)

def worker_extractor(task_queue, result_queue, db_path, ref_fasta_path):
    """
    Consumer worker function. Pulls batches of read alignments from the queue, 
    retrieves their raw signal from POD5 files, extracts the corresponding 
    genomic sequence context (K-mer) from the FASTA file, and writes the output.
    """
    # Open a read-only connection to the SQLite index
    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    except sqlite3.OperationalError:
        conn = sqlite3.connect(db_path)
        
    c = conn.cursor()
    
    # Each worker must open its own FASTA handle to avoid multiprocessing conflicts
    ref_fasta = pysam.FastaFile(ref_fasta_path) if ref_fasta_path else None
    
    proc_name = multiprocessing.current_process().name
    temp_filename = f"temp_worker_{proc_name}_{int(time.time())}.tsv"
    
    import pod5
    
    with open(temp_filename, 'w') as out_f:
        while True:
            batch = task_queue.get()
            if batch is None:
                break
                
            for read_data in batch:
                read_id = read_data['id']
                chrom = read_data['chrom']
                aligned_pairs = read_data['pairs']  # List of tuples: (query_pos, ref_pos)
                
                # 1. Look up POD5 path
                c.execute("SELECT path FROM pod5_map WHERE read_name = ?", (read_id,))
                res = c.fetchone()
                if not res: 
                    continue
                
                # 2. Extract full continuous signal
                try:
                    with pod5.Reader(res[0]) as reader:
                        record = reader.read_one(str(read_id))
                        if not record: 
                            continue
                        full_signal = signal_math.get_picoampere_signal(record)
                except Exception as e:
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