"""
Configuration parameters for PoreSignalExtractor.
This module centralizes default values, constants, and batch sizes used across 
different extraction executors to preserve the exact behavior of the original pipelines.
"""

# ==========================================
# Multiprocessing & Queue Settings
# ==========================================

# Maximum size for the multiprocessing task queues. 
# Keeps memory usage bounded preventing Out-Of-Memory (OOM) errors.
QUEUE_MAX_SIZE = 100

# Number of tasks assigned to a child process before it gets restarted 
# (used in Pool.imap_unordered to free up memory during heavy genomic processing).
MAX_TASKS_PER_CHILD = 10


# ==========================================
# Batching & Chunking Settings
# ==========================================

# Number of reads grouped in a single batch sent to worker queues 
# (used in self-mapping modes like uncalled_self and f5c_self).
READ_BATCH_SIZE = 200

# Number of genomic bases per chunk when splitting the BAM file references 
# (used in moves mode to distribute genome regions across cores).
GENOMIC_CHUNK_SIZE = 1_000_000

# Number of rows processed in a single chunk when reading the eventalign TSV 
# (used in uncalled_ref mode to prevent loading massive tables into RAM).
EVENTALIGN_CHUNK_SIZE = 500_000


# ==========================================
# POD5 I/O & Memory Management
# ==========================================

# Maximum number of POD5 file readers to keep open in memory simultaneously.
# Operates as a Least Recently Used (LRU) cache.
POD5_CACHE_CAPACITY = 16


# ==========================================
# Sequence Context (K-mer) Settings
# ==========================================

# Default length for the genomic K-mer extraction.
# A length of 5 corresponds to: 2 bases upstream + 1 target base + 2 bases downstream.
KMER_LENGTH = 5