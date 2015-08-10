import numpy as np

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
        self.RATIO = False

    def census(self, sequence, track, **kwargs):
        raise NotImplementedError("strategy.census")

class NucleotideCounterStrategy(BaseStrategy):

    def __init__(self, bases):
        if len(bases) == 0:
            bases = ['A', 'C', 'G', 'T', 'N']
        super(NucleotideCounterStrategy, self).__init__(tracks=bases, title="Base Count")

    def census(self, sequence, track, **kwargs):
        #return sequence.count(track)  ## doesn't work on buffers...
        count = 0
        for base in sequence:
            if base == track:
                count += 1
        return count

class MotifCounterStrategy(BaseStrategy): # pragma: no cover

    def __init__(self, motifs, overlap=True):
        import re
        self.modules = { "re": re }
        self.overlap = overlap

        super(MotifCounterStrategy, self).__init__(tracks=motifs, title="Motif Count")

    def census(self, sequence, track, **kwargs):
        if self.overlap:
            return len(self.modules["re"].findall("(?=(%s))" % track, sequence))
        else:
            return len(self.modules["re"].findall(track, sequence))

class PositionCounterStrategy(BaseStrategy): # pragma: no cover

    def __init__(self, tracks=None):
        super(PositionCounterStrategy, self).__init__(title="Count")

    def census(self, sequence, track, **kwargs):
        # Populate the chromosome array with 1 for each position a variant exists
        #TODO Assuming 1-index, perhaps use a kwarg
        count = 0
        for variant_loc in data:
            variant_loc -= kwargs['start']
            if variant_loc > 0:
                try:
                    arr[variant_loc-1] = 1
                    count += 1
                except:
                    break
        return count


class GCRatioStrategy(BaseStrategy):

    def __init__(self, tracks=None):
        super(GCRatioStrategy, self).__init__(title="GC Ratio")

        import re
        self.modules = { "re": re }

        self.RATIO = True
        self.RATIO_OF = None

    def census(self, sequence, track, **kwargs):
        return float(len(self.modules["re"].findall("[GC]", sequence)))/len(sequence)

class ReferenceConsensusStrategy(BaseStrategy): # pragma: no cover

    def __init__(self, tracks=None, polarity=1, reference=None):
        self.POLARITY = polarity
        self.REFERENCE = reference

        title = "Reference Concordance"
        if polarity > 0:
            title = "Reference Matches"
        elif polarity < 0:
            title = "Reference Mismatches"

        super(ReferenceConsensusStrategy, self).__init__(title=title)

    def census(self, sequence, track, **kwargs):
        # Currently only handles global references (ie. not for group/track)
        count = 0
        for location, base in enumerate(data):
            if self.POLARITY > 0:
                if base.upper() == self.REFERENCE[kwargs['chrom']][location].upper():
                    count += 1
            elif self.POLARITY < 0:
                if base.upper() != self.REFERENCE[kwargs['chrom']][location].upper():
                    count += 1
            else:
                if base.upper() == self.REFERENCE[kwargs['chrom']][location].upper():
                    count += 1
                else:
                    count -= 1
        return count
