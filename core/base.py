"""
Core utilities for logging setup and POD5 file indexing.
"""
import os
import sqlite3
import logging
import pod5

logger = logging.getLogger(__name__)

def setup_logger():
    """
    Initializes and configures the global logger.
    """
    logging.basicConfig(
        level=logging.INFO, 
        format='%(asctime)s - %(processName)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger()

def index_pod5_files(pod5_dirs, db_path):
    """
    Creates an SQLite database mapping read IDs to their specific POD5 file paths.
    This enables highly efficient, RAM-friendly random access to reads across 
    thousands of POD5 files during multiprocessing.
    
    Args:
        pod5_dirs (list): List of directory paths containing .pod5 files.
        db_path (str): Path where the SQLite database will be saved/created.
    """
    logger.info(f"Indexing POD5 files from {len(pod5_dirs)} directories into {db_path}...")
    
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS pod5_map
                 (read_name TEXT PRIMARY KEY, path TEXT)''')
    
    file_count = 0
    read_count = 0
    
    for pod5_dir in pod5_dirs:
        for root, _, files in os.walk(pod5_dir):
            for f in files:
                if f.endswith(".pod5"):
                    path = os.path.join(root, f)
                    try:
                        with pod5.Reader(path) as reader:
                            read_ids = reader.read_ids
                            # Batch insertion for performance
                            data = [(str(rid), path) for rid in read_ids]
                            c.executemany("INSERT OR REPLACE INTO pod5_map VALUES (?, ?)", data)
                            read_count += len(data)
                        file_count += 1
                    except Exception as e:
                        logger.warning(f"Error reading POD5 file {path}: {e}")
    
    conn.commit()
    
    # Enable Write-Ahead Logging (WAL) for better concurrent read performance across workers
    c.execute("PRAGMA journal_mode=WAL") 
    conn.close()
    
    logger.info(f"Indexing complete: {read_count} reads mapped across {file_count} POD5 files.")