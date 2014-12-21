import numpy as np

class StrategyValue(float):

    def __new__(cls, value, *args, **kwargs):
        return float.__new__(cls, value)

    def __init__(self, value, k=None):
        # Initialise default weight if one has not been provided by the user
        # Any non-zero value is weighted 1 otherwise it is 0 so as to prevent
        # StrategyValue's initialised to 0 from having any weighting.
        if k is None:
            if value == 0.0:
                k = 0.0
            else:
                k = 1.0

        self.k = k

    def __add__(self, other):
        #NOTE Catch AttributeError for 'other' objects that have no k
        #     In these cases, new_k is always 1, i.e. there is no weighting
        try:
            current_total = self * self.k
            other_total = other * other.k

            new_k = self.k + other.k
            new_average = (current_total + other_total) / new_k
        except AttributeError:
            new_k = 1
            new_average = (float(self) + float(other)) / new_k

        return StrategyValue(new_average, new_k)

class BaseStrategy(object):

    def __init__(self, tracks=None, title=""):
        if tracks is None:
            tracks = ["default"]

        self.TRACKS = tracks
        self.AXIS_TITLE = title

    def prepare(self, arr, data, track):
        """Populate elements in a given iterable 'arr' (typically a numpy
        array) following some processing strategy on genomic sequence 'data'.
        Tracks can be used to further inform the strategy on what
        behaviour to use (for example a nucleotide counting strategy would need
        a 'track' for each base to be counted)."""
        raise NotImplementedError("strategy.prepare")

    def evaluate(self, region, **kwargs):
        """Evaluate the contents of a given iterable 'region' (typically a numpy
        array) as prepared by this strategy. The simplest strategies will sum
        the binary flags over the array."""
        raise NotImplementedError("strategy.evaluate")

class NucleotideCounterStrategy(BaseStrategy):

    def __init__(self, bases):
        super(NucleotideCounterStrategy, self).__init__(tracks=bases, title="Base Count")

    def prepare(self, arr, data, current_track):
        for location, base in enumerate(data):
            if base.upper() == current_track:
                arr[location] = 1
        return arr

    def evaluate(self, region, **kwargs):
        return np.sum(region)

class KMerCounterStrategy(BaseStrategy): # pragma: no cover

    def __init__(self, kmers):
        super(KMerCounterStrategy, self).__init__(tracks=kmers, title="Motif Count")

    def prepare(self, arr, data, track):
        import re

        # Populate the region array with 1 for the start position of desired K-Mer
        #TODO Use KMP?
        for location in [m.start() for m in re.finditer(track, data)]:
            arr[location] = 1
        return arr

    def evaluate(self, region, **kwargs):
        # Need to allow for cases where the K-mer present flag has been placed
        # at the very end of the region (and thus the region merely contains only
        # the first base of the desired K-mer
        return np.sum(region[:-(len(kwargs['track'])-1)])

class VariantCounterStrategy(BaseStrategy): # pragma: no cover

    def __init__(self, tracks=None):
        super(VariantCounterStrategy, self).__init__(title="Variant Count")

    def prepare(self, arr, data, track):
        """Return a NumPy array containing 1 for position elements where a variant
        exists and 0 otherwise."""
        # Populate the chromosome array with 1 for each position a variant exists
        for variant_loc in data:
            arr[variant_loc] = 1
        return arr

    def evaluate(self, region, **kwargs):
        return np.sum(region)


class GCRatioStrategy(BaseStrategy): # pragma: no cover

    def __init__(self, tracks=None):
        super(GCRatioStrategy, self).__init__(title="GC Ratio")

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

    def evaluate(self, region, **kwargs):
        return StrategyValue(float(np.sum(region))/len(region), len(region))
