"""
Executor for extracting raw signal based on Uncalled4 reference event alignment.
Processes large TSV mapping files using Pandas chunking, distributes chunks to 
workers via multiprocessing, and uses optimized batch I/O for POD5 extraction.
"""
import os
import time
import pandas as pd
import sqlite3
import shutil
import multiprocessing
import logging

from core import config
from core import base
from core import signal_math
from core import bed_handling

logger = logging.getLogger(__name__)

def worker_extractor(task_queue, result_queue, db_path):
    """
    Consumer worker function. Pulls Pandas dataframes (chunks of the TSV) from the queue, 
    groups the required reads by POD5 file, fetches the signal in batch, and writes the output.
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
            chunk = task_queue.get()
            if chunk is None:
                break
                
            if chunk.empty:
                continue

            # 1. Trova le read uniche nel chunk
            unique_reads = chunk['read_name'].unique()
            
            # 2. Raggruppa le read per file POD5 (Ottimizzazione batch originale)
            file_batch_map = {}
            for r_id in unique_reads:
                str_r_id = str(r_id)
                c.execute("SELECT path FROM pod5_map WHERE read_name = ?", (str_r_id,))
                res = c.fetchone()
                if res:
                    p_path = res[0]
                    if p_path not in file_batch_map:
                        file_batch_map[p_path] = []
                    file_batch_map[p_path].append(str_r_id)

            # 3. Estrae i segnali in blocco (batch reading) in RAM
            signal_cache = {}
            for p_path, read_ids in file_batch_map.items():
                try:
                    with pod5.Reader(p_path) as reader:
                        for record in reader.reads(read_ids):
                            signal_cache[str(record.read_id)] = signal_math.get_picoampere_signal(record)
                except Exception as e:
                    continue

            # 4. Processa le righe del chunk di Pandas ed estrae le fette
            grouped = chunk.groupby("read_name")
            buffer = []
            
            for read_id, group in grouped:
                read_id_str = str(read_id)
                
                if read_id_str not in signal_cache:
                    continue
                
                full_signal = signal_cache[read_id_str]
                sig_len = len(full_signal)
                
                for row in group.itertuples():
                    start = int(row.start_idx)
                    end = int(row.end_idx)
                    
                    if start < 0 or end > sig_len or start >= end:
                        continue
                    
                    raw_slice = full_signal[start:end]
                    sig_str = ",".join([f"{x:.1f}" for x in raw_slice])
                    
                    # Usa il formato a 6 colonne (inserendo '.' per read_kmer)
                    buffer.append(f"{read_id}\t{row.contig}\t{row.position}\t{row.reference_kmer}\t.\t{sig_str}")
            
            if buffer:
                out_f.write("\n".join(buffer) + "\n")
                
            # Libera la RAM
            del signal_cache
            del file_batch_map

    conn.close()
    result_queue.put(temp_filename)


def run(args):
    """
    Main entry point for the 'unc-ref' extraction subcommand.
    """
    db_path = "uncalled_ref_map.db"
    
    # 0. Load BED target sites
    bed_target_sites = bed_handling.load_bed_target_sites(getattr(args, 'bed', None))
    
    # 1. Indexing
    logger.info("Indexing POD5 files into SQLite database...")
    if os.path.exists(db_path):
        os.remove(db_path)
    base.index_pod5_files(args.pod5, db_path)
    
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
            args=(task_queue, result_queue, db_path)
        )
        p.start()
        workers.append(p)
    
    # 4. Producer: Read TSV mapping file in chunks
    logger.info(f"Parsing mapping file: {args.map} in chunks of {config.EVENTALIGN_CHUNK_SIZE} rows")
    
    map_iter = pd.read_csv(
        args.map, 
        sep='\t', 
        chunksize=config.EVENTALIGN_CHUNK_SIZE, 
        usecols=['read_name', 'contig', 'position', 'reference_kmer', 'start_idx', 'end_idx']
    )
    
    chunk_count = 0
    for chunk in map_iter:
        # --- APPLICAZIONE DEL FILTRO BED NEL PRODUCER ---
        # Filtriamo il chunk prima di inviarlo ai worker per risparmiare tempo di comunicazione
        if bed_target_sites is not None:
            mask = chunk.apply(lambda r: bed_target_sites.contains(r.contig, r.position), axis=1)
            chunk = chunk[mask]
        
        if not chunk.empty:
            task_queue.put(chunk)
            
        chunk_count += 1
        if chunk_count % 10 == 0:
            logger.info(f"Dispatched {chunk_count} chunks to workers...")
            
    # Send poison pills
    for _ in range(args.threads):
        task_queue.put(None)
        
    logger.info("TSV parsing completed. Waiting for workers to finish...")

    # 5. Collect Results
    temp_files = [result_queue.get() for _ in range(args.threads)]
    
    for p in workers:
        p.join()

    # 6. Merge Files
    logger.info("Merging temporary worker files into final output...")
    with open(args.output, 'w') as out_f:
        out_f.write("read_id\tcontig\tposition\tref_kmer\tread_kmer\tsamples_pA\n")
        
        for tf in temp_files:
            if os.path.exists(tf):
                with open(tf, 'r') as in_f:
                    shutil.copyfileobj(in_f, out_f)
                os.remove(tf)

    # Cleanup SQLite DB
    if os.path.exists(db_path):
        os.remove(db_path)
        
    logger.info("Uncalled4 reference extraction completed successfully.")
