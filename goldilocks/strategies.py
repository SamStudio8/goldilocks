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

