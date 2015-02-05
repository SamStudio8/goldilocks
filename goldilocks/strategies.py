import numpy as np

class StrategyValue(float):

    def __new__(cls, value, *args, **kwargs):
        return float.__new__(cls, value)

    def __init__(self, value, k=0):
        if k != 0:
            if k < 1:
                raise ValueError("k must be == 0 or >= 1")
        self.k = k

    def __add__(self, other):

        # Try to check weight of other, floats and the like will
        # raise an AttributeError
        try:
            other_k = other.k
        except AttributeError:
            other_total = float(other)
            other_k = 0

        # Weight other
        if other_k >= 1:
            other_total = other * other.k
        else:
            other_total = float(other)

        # Weight self
        if self.k >= 1:
            current_total = self * self.k
        else:
            current_total = float(self)

        # New weight k
        new_k = self.k + other_k
        if new_k == 0:
            new_average = (current_total + other_total)
        else:
            new_average = (current_total + other_total) / new_k

        return StrategyValue(new_average, new_k)

    def __radd__(self, other):
        return self.__add__(other)

class BaseStrategy(object):
    """Interface for which census strategies must be compliant.

    It is intended that all valid census strategies must inherit from `BaseStrategy`
    and provide implementations for each of the methods.

    Parameters
    ----------
    tracks : list{str}, optional(default=None)
        A list of strings defining multiple features of interest in the context
        of this strategy with which to perform the census. For example a simple
        nucleotide counting strategy will accept a list of nucleotides of interest.

        By default the argument is None which will cause the `TRACKS` attribute
        to be populated with one track; "default".

    title : str, optional(default="")
        A string used to annotate what the values returned by this strategy
        `evaluate` method represent, particularly for use on the y-axis in any
        plots generated.

        If performing a census for GC Ratio, a suitable title might be "GC Ratio".
        A nucleotide counter may use something more generic such as "Nucleotide Count".

    Attributes
    ----------
    TRACKS : list{str}
        A list of strings representing features of interest in the input data
        to be censused as provided by the user on instantiation.

    AXIS_TITLE : str
        A string used to establish the context of the values returned from this
        strategy's `evaluate` function, as provided by the user on instantiation.
    """
    def __init__(self, tracks=None, title=""):
        if tracks is None:
            tracks = ["default"]

        self.TRACKS = tracks
        self.AXIS_TITLE = title

    def prepare(self, arr, data, track):
        """Parse genomic data and apply some algorithm ('strategy') to populate
        an array with values for later evaluation.

        Populate elements in a given iterable `arr` (typically a numpy array)
        by applying some processing strategy on genomic sequence `data`.
        `track` can be used to further inform this function's strategy on what
        behaviour to apply. For example, a nucleotide counting strategy would
        use `track` to discern which nucleotide to count in `data`.

        Parameters
        ----------
        arr : array_like

        data : array_like

        track : str

        Returns
        -------
        arr : array_like

        Raises
        ------
        NotImplementedError
            If attempting to utilise a strategy that inherits from `BaseStrategy`
            without providing an implementation for `prepare`.
        """
        raise NotImplementedError("strategy.prepare")

    def evaluate(self, region, **kwargs):
        """Evaluate the contents of a given iterable 'region' (typically a numpy
        array) as prepared by this strategy. The simplest strategies will sum
        the binary flags over the array.

        Raises
        ------
        NotImplementedError
            If attempting to utilise a strategy that inherits from `BaseStrategy`
            without providing an implementation for `evaluate`.
        """
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
        # Populate the chromosome array with 1 for each position a variant exists
        for variant_loc in data:
            arr[variant_loc] = 1
        return arr

    def evaluate(self, region, **kwargs):
        return np.sum(region)


class GCRatioStrategy(BaseStrategy):

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
