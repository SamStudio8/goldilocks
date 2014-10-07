import numpy as np

class KMerCounterStrategy(object):

    def __init__(self, kmers):
        self.TRACKS = kmers

    def prepare(self, arr, data, track):
        import re

        # Populate the region array with 1 for the start position of desired K-Mer
        #TODO Use KMP?
        for location in [m.start() for m in re.finditer(track, data)]:
            arr[location] = 1
        return arr

    def evaluate(self, region, track):
        # Need to allow for cases where the K-mer present flag has been placed
        # at the very end of the region (and thus the region merely contains only
        # the first base of the desired K-mer
        return np.sum(region[:-(len(track)-1)])

class VariantCounterStrategy(object):

    def __init__(self, tracks=None):
        self.TRACKS = ["default"]

    def prepare(self, arr, data, track):
        """Return a NumPy array containing 1 for position elements where a variant
        exists and 0 otherwise."""
        # Populate the chromosome array with 1 for each position a variant exists
        for variant_loc in data:
            arr[variant_loc] = 1
        return arr

    def evaluate(self, region, track):
        return np.sum(region)


class GCRatioStrategy(object):

    def __init__(self, tracks=None):
        self.TRACKS = ["default"]

    def prepare(self, arr, data, track):
        # Populate the region array with 1 for each position a GC base exists
        #for location, base in enumerate(data):
        #    base_u = base.upper()
        #    if base_u == "G" or base_u == "C":
        #        arr[location] = 1
        #return arr
        import re
        for location in [m.start() for m in re.finditer("[GC]", data)]:
            arr[location] = 1
        return arr

    def evaluate(self, region, track):
        return float(np.sum(region))/len(region)


class NCounterStrategy(object):

    def __init__(self, tracks=None):
        self.TRACKS = ["default"]

    def prepare(self, arr, data, track):
        # Populate the region array with 1 for each position a missing base exists
        for location, base in enumerate(data):
            if base.upper() == "N":
                arr[location] = 1
        return arr

    def evaluate(self, region, track):
        return np.sum(region)

