=================
Custom Strategies
=================

One of the major features of Goldilocks is its extensibility. Strategies
are both easily customisable and interchangeable, as they all share a
common interface. This interface also provides a platform for users with
some knowledge of Python to construct their own custom census rules. One
such example follows below:

A Simple ORF Finder
-------------------

Code Sample
~~~~~~~~~~~

::

    # Import Goldilocks and the BaseStrategy class
    from goldilocks import Goldilocks
    from goldilocks.strategies import BaseStrategy

    # Define a new class for your custom strategy that inherits from BaseStrategy
    class MyCustomSimpleORFCounterStrategy(BaseStrategy):
        
        # Initialising function boilerplate, required to set-up some properties of the census
        def __init__(self, tracks=None, min_codons=1):
            # Initialise the custom class with super
            super(MyCustomSimpleORFCounterStrategy, self).__init__(
                    tracks=range(0,3),                   # Use range to specify a counter for
                                                        #  each of the three possible forward
                                                        #  reading frames in which to search
                                                        #  to search for open reading frames
                    label="Forward Open Reading Frames"  # Y-Axis Plot Label
            )
            self.MIN_CODONS = min_codons

        # This function defines the actual behaviour of a census for a given region
        #  of sequence and the current counting track (one of three reading frames)
        def census(self, sequence, track_frame, **kwargs):
            STARTS = ["ATG"]
            STOPS = ["TAA", "TGA", "TAG"]
            CODON_SIZE = 3

            # Split input sequence into codons. Open a frame if a START is found
            #  and increment the ORF counter if a STOP is encountered afterward
            orfs = orf_open = 0
            for i in xrange(track_frame, len(sequence), CODON_SIZE):
                codon = sequence[i:i+CODON_SIZE].upper()
                if codon in STARTS and orf_open == 0:
                    orf_open = 1
                elif codon in STOPS and orf_open > 0:
                    if orf_open > self.MIN_CODONS:
                        orfs += 1
                    orf_open = 0
                elif orf_open > 0:
                    orf_open += 1
            return orfs

    # Organise and execute the census
    sequence_data = { "hs37d5": {"file": "/store/ref/hs37d5.1-3.fa.fai"} }
    g = Goldilocks(MyCustomSimpleORFCounterStrategy(min_codons=30), sequence_data,
            length="1M", stride="1M", is_faidx=True, processes=4)

Implementation Description
~~~~~~~~~~~~~~~~~~~~~~~~~~

Strategies are defined as Python classes, inheriting from the
``BaseStrategy`` class found in the ``goldilocks.strategies``
subpackage. The class requires just two function definitions to be
compliant with the shared interface; ``__init__``: the class initializer
that takes care of the setup of the strategy’s internals via the
``BaseStrategy`` parent class, and ``census``: the function actually
responsible for the behaviour of the strategy itself.

The example presented is a very simple open reading frame counter. It
searches the three forward frames for start codons that are then
followed by one of the three stop codons. The \`\`tracks" in this
example are the three possible frames. Note on line 9 that our
``__init__`` provides a default argument for ``tracks`` of ``None``.
Thus this particular strategy does not need the ``tracks`` argument.
Instead, the track list is provided by the strategy itself, and passed
to the ``BaseStrategy`` ``__init__`` (line 12), forcing tracks to be the
list [0, 1, 2]. The elements of this list are used as an integer offset
from which to begin splitting input DNA sequences when conducting the
census later, which is why on this occasion we don’t want to allow the
user to specify their own tracks. Other strategies, such as the included
``NucleotideCounterStrategy`` just pass the ``tracks`` argument from the
user through to the super ``__init__``.

For a given array of ``sequence`` data and a frame offset
(``track_frame``), the ``census`` function splits the sequence into
nucleotide triplets from the offset and searches for open reading
frames. A subsequence is considered an ORF by this strategy if the ATG
START codon is encountered and later followed by any STOP codon.

Our example finishes with the familiar specification of the location of
input sequence data and the construction of the census itself. Here we
specify a census of all 1Mbp regions with no overlap (that is, the
stride is equal to the size of the regions) and instantiate our new
``MyCustomSimpleORFCounterStrategy`` with a keyword requiring valid ORFs
to be at least 30 codons in length (excluding start and stop).

Every strategy’s ``census`` function is expected to return a numerical
result that can be used to rank and sort regions, in this scenario,
``census`` returns the number of ORFs found.

Note also, strategies may specify any number of keyword arguments that
are not found in the ``BaseStrategy``. In our example, ``min_codons``
can be set by a user to specify how many codons must lie between an
opening and closing codon to be counted as an open reading frame. We
store this value as a member of the strategy object on line 18 and use
it on line 35 to ensure the ``orfs`` counter is only incremented when
the length of the current open reading frame has exceeded the provided
threshold. One could store any number of configurable parameters inside
of the strategy class in this fashion. This framework allows one to
increase the complexity of strategies while still providing a friendly
and interchangeable interface for end users.
