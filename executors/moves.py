"""
Executor for extracting raw signal using Dorado's implicit alignment (BAM moves).
This module parallelizes the extraction by splitting the reference genome into chunks,
processing multiple regions concurrently while managing POD5 memory footprints via an LRU cache.
"""
import os
import pysam
import numpy as np
import logging
from multiprocessing import Pool
from functools import partial
import shutil
import gc
import sqlite3

from core import config
from core import signal_math
from core import base
from core import file_management
from core import bed_handling

logger = logging.getLogger(__name__)

# --- GLOBAL WORKER STATE ---
worker_cache = None
worker_db_conn = None
worker_fasta = None
worker_bed_sites = None

def init_worker(db_path, ref_file, bed_sites):
    """
    Initializes resources for each worker process in the multiprocessing pool.
    Sets up the POD5 file cache, the SQLite connection, FASTA reader, and BED targets.
    """
    global worker_cache, worker_db_conn, worker_fasta, worker_bed_sites
    
    worker_cache = file_management.Pod5FileCache()
    worker_bed_sites = bed_sites
    
    # Open SQLite in read-only mode for thread-safe concurrent access
    db_uri = f"file:{db_path}?mode=ro"
    worker_db_conn = sqlite3.connect(db_uri, uri=True)
    
    if ref_file:
        worker_fasta = pysam.FastaFile(ref_file)

def process_region(region_task, bam_path, output_dir, include_unmapped):
    """
    Processes a specific genomic region (or unmapped reads), extracting the raw signal 
    segments based on the basecaller's 'mv' (moves) tags.
    """
    global worker_cache, worker_db_conn, worker_fasta, worker_bed_sites
    
    chrom, start_pos, end_pos, task_id = region_task
    temp_filename = os.path.join(output_dir, f"temp_{task_id}.tsv")
    c = worker_db_conn.cursor()
    
    try:
        bam = pysam.AlignmentFile(bam_path, "rb", check_sq=False)
    except Exception as e:
        logger.error(f"Task {task_id}: Error opening BAM {bam_path} -> {e}")
        return None

    if chrom == "UNMAPPED":
        try:
            iterator = bam.fetch(until_eof=True)
            iterator = (r for r in iterator if r.is_unmapped)
        except Exception as e:
            return None
    else:
        try:
            iterator = bam.fetch(chrom, start_pos, end_pos)
        except Exception as e:
            return None

    reads_written = 0
    kmer_offset = config.KMER_LENGTH // 2

    with open(temp_filename, "w") as out_f:
        for read in iterator:
            read_id = read.query_name

            if not read.has_tag('mv'): 
                continue
            
            # 1. Get POD5 path
            c.execute("SELECT path FROM pod5_map WHERE read_name = ?", (read_id,))
            res = c.fetchone()
            if not res: 
                continue
            
            pod5_path = res[0]

            try:
                moves_info = read.get_tag('mv')
                stride = moves_info[0]
                moves_array = np.array(moves_info[1:], dtype=np.uint8)
                
                # 2. Get reader from LRU cache
                reader = worker_cache.get_reader(pod5_path)
                if not reader: 
                    continue
                
                try:
                    record = reader.read_one(read_id)
                    if record is None: 
                        continue
                    full_signal = signal_math.get_picoampere_signal(record)
                except Exception:
                    continue

                ts_shift = read.get_tag('ts') if read.has_tag('ts') else 0
                seq_len = read.query_length
                
                valid_moves_indices = np.where(moves_array == 1)[0]
                
                if len(valid_moves_indices) < seq_len:
                    continue

                base_to_move_idx = valid_moves_indices[:seq_len]
                
                if not read.is_unmapped:
                    aligned_pairs = read.get_aligned_pairs(matches_only=False)
                    read_to_ref_map = {q: r for q, r in aligned_pairs if q is not None}
                else:
                    read_to_ref_map = {}

                read_seq = read.query_sequence
                chrom_out = read.reference_name if not read.is_unmapped else "UNMAPPED"

                # Iterate through bases
                for base_idx in range(kmer_offset, seq_len - kmer_offset):
                    ref_pos = read_to_ref_map.get(base_idx)
                    
                    # --- FILTRO BED ---
                    if worker_bed_sites is not None:
                        # Se è unmapped o la posizione non è nel BED, salta questa base
                        if ref_pos is None or not worker_bed_sites.contains(chrom_out, ref_pos):
                            continue
                    # ------------------
                    
                    start_move_idx = base_to_move_idx[base_idx]
                    end_move_idx = base_to_move_idx[base_idx + 1] if base_idx + 1 < seq_len else len(moves_array)
                    
                    start_signal_idx = (start_move_idx * stride) + ts_shift
                    end_signal_idx = (end_move_idx * stride) + ts_shift

                    if start_signal_idx >= len(full_signal) or start_signal_idx >= end_signal_idx:
                        continue
                        
                    segment_signal = full_signal[start_signal_idx:end_signal_idx]
                    read_kmer = read_seq[base_idx - kmer_offset : base_idx + kmer_offset + 1]
                    
                    if ref_pos is not None and worker_fasta is not None:
                        ref_kmer = signal_math.get_kmer_from_fasta(worker_fasta, chrom_out, ref_pos)
                        pos_out = ref_pos
                    else:
                        ref_kmer = "N" * config.KMER_LENGTH
                        pos_out = "NA"

                    signal_str = ",".join(f"{x:.1f}" for x in segment_signal)
                    
                    # Formato a 6 colonne standard
                    row_data = [
                        read_id, 
                        chrom_out, 
                        str(pos_out), 
                        ref_kmer, 
                        read_kmer,
                        signal_str
                    ]
                        
                    out_f.write("\t".join(row_data) + "\n")
                    reads_written += 1
                    
            except Exception as e:
                continue

    bam.close()
    
    if task_id % 50 == 0:
        gc.collect()

    return temp_filename

def run(args):
    """
    Main entry point for the 'moves' extraction subcommand.
    """
    db_path = "moves_map.db"
    
    # 0. Load BED target sites
    bed_target_sites = bed_handling.load_bed_target_sites(getattr(args, 'bed', None))
    
    # 1. Indexing
    base.index_pod5_files(args.pod5, db_path)

    out_dir = "temp_signal_chunks"
    os.makedirs(out_dir, exist_ok=True)

    logger.info("Reading BAM header to generate genomic tasks...")
    bam = pysam.AlignmentFile(args.bam, "rb", check_sq=False)
    ref_lens = {chrom: bam.get_reference_length(chrom) for chrom in bam.references}
    bam.close()

    # 2. Task Generation (Chuncking del genoma per parallelizzazione su interi cromosomi)
    tasks = []
    pid_counter = 0
    
    if args.include_unmapped:
        tasks.append(("UNMAPPED", 0, 0, pid_counter))
        pid_counter += 1
        
    if ref_lens:
        for chrom, length in ref_lens.items():
            for start in range(0, length, config.GENOMIC_CHUNK_SIZE):
                end = min(start + config.GENOMIC_CHUNK_SIZE, length)
                tasks.append((chrom, start, end, pid_counter))
                pid_counter += 1

    logger.info(f"Generated {len(tasks)} tasks. Starting Multiprocessing Pool ({args.threads} cores).")
    
    func = partial(
        process_region,
        bam_path=args.bam,
        output_dir=out_dir,
        include_unmapped=args.include_unmapped
    )

    # 3. Execution
    temp_files = []
    with Pool(
        processes=args.threads, 
        initializer=init_worker, 
        initargs=(db_path, getattr(args, 'fasta', None), bed_target_sites), 
        maxtasksperchild=config.MAX_TASKS_PER_CHILD
    ) as pool:
        for res in pool.imap_unordered(func, tasks, chunksize=1):
            if res: 
                temp_files.append(res)

    # 4. Merging
    logger.info("Merging temporary files into final output...")
    # Header a 6 colonne standardizzato!
    header = ['read_id', 'contig', 'position', 'ref_kmer', 'read_kmer', 'samples_pA']
    
    with open(args.output, 'w') as final_out:
        final_out.write('\t'.join(header) + '\n')
        for f in temp_files:
            if os.path.exists(f):
                with open(f, 'r') as temp_in:
                    shutil.copyfileobj(temp_in, final_out)
                os.remove(f)
                
    shutil.rmtree(out_dir)
    if os.path.exists(db_path):
        os.remove(db_path)
        
    logger.info("Moves extraction pipeline completed successfully.")
