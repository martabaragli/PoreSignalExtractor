import logging
import bisect

logger = logging.getLogger(__name__)

class BedRegions:
    def __init__(self, bed_path):
        self.regions = {}
        self._load_and_merge(bed_path)

    def _load_and_merge(self, bed_path):
        logger.info(f"Loading target regions from BED file: {bed_path}")
        raw_regions = {}
        
        with open(bed_path) as f:
            for line in f:
                if line.startswith("#") or not line.strip(): 
                    continue
                parts = line.strip().split()
                if len(parts) >= 3:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    if chrom not in raw_regions:
                        raw_regions[chrom] = []
                    raw_regions[chrom].append((start, end))
        
        # Fonde gli intervalli sovrapposti per una ricerca iper-veloce
        for chrom, intervals in raw_regions.items():
            intervals.sort(key=lambda x: x[0])
            merged = []
            for interval in intervals:
                if not merged or merged[-1][1] < interval[0]:
                    merged.append(interval)
                else:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], interval[1]))
            
            # Salva liste separate per start e end (ottimizzate per bisect)
            self.regions[chrom] = {
                'starts': [m[0] for m in merged],
                'ends': [m[1] for m in merged]
            }

    def contains(self, chrom, pos):
        """Restituisce True se la posizione cade in una delle regioni per il cromosoma dato."""
        if chrom not in self.regions:
            return False
            
        starts = self.regions[chrom]['starts']
        ends = self.regions[chrom]['ends']
        
        # Ricerca binaria O(log N) - Estremamente veloce
        idx = bisect.bisect_right(starts, pos) - 1
        if idx >= 0 and pos < ends[idx]:
            return True
        return False

def load_bed_target_sites(bed_path):
    """Restituisce l'oggetto BedRegions o None se non è fornito un file BED."""
    if not bed_path:
        return None
    return BedRegions(bed_path)
