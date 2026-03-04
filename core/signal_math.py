"""
Mathematical and sequence-related utility functions for signal processing.
"""
import math
import logging
from core import config

logger = logging.getLogger(__name__)

def get_picoampere_signal(read_record):
    """
    Converts the raw integer signal array from a POD5 record 
    into calibrated picoampere (pA) measurements.
    
    Args:
        read_record (pod5.ReadRecord): A single read record from a POD5 file.
        
    Returns:
        numpy.ndarray: Array of calibrated pA signal values.
    """
    cal = read_record.calibration
    return (read_record.signal + cal.offset) * cal.scale


def get_kmer_from_fasta(fasta_handle, chrom, pos, k=config.KMER_LENGTH):
    """
    Fetches the sequence context (k-mer) around a given genomic position.
    
    Args:
        fasta_handle (pysam.FastaFile): Open handle to the reference FASTA.
        chrom (str): Chromosome name.
        pos (int): Genomic position (0-based).
        k (int): K-mer length (must be odd). Defaults to config.KMER_LENGTH.
        
    Returns:
        str: Extracted k-mer sequence in uppercase, or Ns if fetching fails.
    """
    try:
        offset = math.floor(k / 2)
        start = pos - offset
        end = pos + offset + 1
        kmer = fasta_handle.fetch(chrom, start, end)
        return kmer.upper() if kmer else "N" * k
    except Exception as e:
        logger.debug(f"Failed to fetch k-mer at {chrom}:{pos} -> {e}")
        return "N" * k