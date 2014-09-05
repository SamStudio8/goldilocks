import numpy as np

class KMerCounterStrategy(object):

    def __init__(self, kmers):
        self.TRACKS = kmers

    def prepare(self, size, data, track):
        import re

        # Populate the region array with 1 for the start position of desired K-Mer
        chro = np.zeros(size+1, np.int8)

        #TODO Use KMP?
        for location in [m.start() for m in re.finditer(track, data)]:
            chro[location] = 1
        return chro

    def evaluate(self, region, track):
        # Need to allow for cases where the K-mer present flag has been placed
        # at the very end of the region (and thus the region merely contains only
        # the first base of the desired K-mer
        return np.sum(region[:-(len(track)-1)])

class VariantCounterStrategy(object):

    def __init__(self, tracks=None):
        self.TRACKS = ["1"]

    def prepare(self, size, data, track):
        """Return a NumPy array containing 1 for position elements where a variant
        exists and 0 otherwise."""
        chro = np.zeros(size+1, np.int8)

        # Populate the chromosome array with 1 for each position a variant exists
        for variant_loc in data:
            chro[variant_loc] = 1
        return chro

    def evaluate(self, region, track):
        return np.sum(region)


class GCRatioStrategy(object):

    def __init__(self, tracks=None):
        self.TRACKS = ["1"]

    def prepare(self, size, data, track):
        # Populate the region array with 1 for each position a GC base exists
        chro = np.zeros(size+1, np.int8)
        for location, base in enumerate(data):
            base_u = base.upper()
            if base_u == "G" or base_u == "C":
                chro[location] = 1
        return chro

    def evaluate(self, region, track):
        return float(np.sum(region))/len(region)


class NCounterStrategy(object):

    def __init__(self, tracks=None):
        self.TRACKS = ["1"]

    def prepare(self, size, data, track):
        # Populate the region array with 1 for each position a missing base exists
        chro = np.zeros(size+1, np.int8)
        for location, base in enumerate(data):
            if base.upper() == "N":
                chro[location] = 1
        return chro

    def evaluate(self, region, track):
        return np.sum(region)

