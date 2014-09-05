import numpy as np

class KMerCounterStrategy(object):

    def __init__(self, kmer):
        self.KMER = kmer

    def prepare(self, size, data):
        import re

        # Populate the region array with 1 for the start position of desired K-Mer
        chro = np.zeros(size+1, np.int8)

        #TODO Use KMP?
        for location in [m.start() for m in re.finditer(self.KMER, data)]:
            chro[location] = 1
        return chro

    def evaluate(self, region):
        # Need to allow for cases where the K-mer present flag has been placed
        # at the very end of the region (and thus the region merely contains only
        # the first base of the desired K-mer
        return np.sum(region[:-(len(self.KMER)-1)])

class VariantCounterStrategy(object):

    @staticmethod
    def prepare(size, data):
        """Return a NumPy array containing 1 for position elements where a variant
        exists and 0 otherwise."""
        chro = np.zeros(size+1, np.int8)

        # Populate the chromosome array with 1 for each position a variant exists
        for variant_loc in data:
            chro[variant_loc] = 1
        return chro

    @staticmethod
    def evaluate(region):
        return np.sum(region)


class GCRatioStrategy(object):

    @staticmethod
    def prepare(size, data):
        # Populate the region array with 1 for each position a GC base exists
        chro = np.zeros(size+1, np.int8)
        for location, base in enumerate(data):
            base_u = base.upper()
            if base_u == "G" or base_u == "C":
                chro[location] = 1
        return chro

    @staticmethod
    def evaluate(region):
        return float(np.sum(region))/len(region)


class NCounterStrategy(object):

    @staticmethod
    def prepare(size, data):
        # Populate the region array with 1 for each position a missing base exists
        chro = np.zeros(size+1, np.int8)
        for location, base in enumerate(data):
            if base.upper() == "N":
                chro[location] = 1
        return chro

    @staticmethod
    def evaluate(region):
        return np.sum(region)

