import numpy as np

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
            if base.upper() == "G" or base.upper() == "C":
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

