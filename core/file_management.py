"""
File management utilities, specifically optimized for handling massive numbers 
of POD5 files without exhausting system memory or file descriptors.
"""
import logging
from collections import deque
import pod5
from core import config

logger = logging.getLogger(__name__)

class Pod5FileCache:
    """
    A Least Recently Used (LRU) cache for POD5 file readers.
    Opening and closing POD5 files repeatedly is computationally expensive.
    This class maintains a pool of active readers up to a specified capacity.
    When the capacity is reached, the least recently accessed file is safely closed.
    """
    def __init__(self, capacity=config.POD5_CACHE_CAPACITY):
        self.capacity = capacity
        self.cache = {} 
        self.queue = deque() 

    def get_reader(self, path):
        """
        Retrieves an active pod5.Reader for the specified path.
        If the reader is already open, it updates its position in the LRU queue.
        If not, it opens a new reader, potentially evicting the oldest one if at capacity.
        
        Args:
            path (str): The absolute or relative path to the POD5 file.
            
        Returns:
            pod5.Reader: An open reader object, or None if the file could not be opened.
        """
        # If the file is already open, mark it as the most recently used
        if path in self.cache:
            self.queue.remove(path)
            self.queue.append(path)
            return self.cache[path]
        
        # If capacity is reached, close and remove the oldest reader
        if len(self.cache) >= self.capacity:
            oldest_path = self.queue.popleft()
            try:
                self.cache[oldest_path].close()
            except Exception as e:
                logger.debug(f"Error closing old POD5 reader {oldest_path}: {e}")
            del self.cache[oldest_path]
        
        # Open the new reader, add it to the cache and the queue
        try:
            reader = pod5.Reader(path)
            self.cache[path] = reader
            self.queue.append(path)
            return reader
        except Exception as e:
            logger.warning(f"Unable to open POD5 file at {path}: {e}")
            return None

    def close_all(self):
        """
        Safely closes all currently open POD5 readers and clears the cache.
        It is crucial to call this before shutting down workers to prevent memory leaks.
        """
        for path, reader in self.cache.items():
            try:
                reader.close()
            except Exception as e:
                logger.debug(f"Error closing POD5 reader {path}: {e}")
                
        self.cache.clear()
        self.queue.clear()