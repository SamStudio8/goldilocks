===========
Basic Usage
===========

The following example assumes basic Python programming experience (and
that you have installed Goldilocks), skip to the
end if you think you know what you're doing.

Importing
---------

To use Goldilocks you will need to import the :class:`goldilocks.goldilocks.Goldilocks`
class and your desired census strategy (e.g. NucleotideCounterStrategy) from
:mod:`goldilocks.strategies` to your script: ::

    from goldilocks.goldilocks import Goldilocks
    from goldilocks.strategies import NucleotideCounterStrategy


Providing Sequence Data
-----------------------

The :class:`goldilocks.goldilocks.Goldilocks` class expects you to provide
sequence data in the following format: ::

    sequence_data = {
        "sample_name_or_identifier": {
            "chr_name_or_number": "my_actual_sequence",
        }
    }

For example: ::

    sequence_data = {
        "my_sample": {
            2: "NANANANANA",
            "one": "CATCANCAT",
            "X": "GATTACAGATTACAN"
        },
        "my_other_sample": {
            2: "GANGANGAN",
            "one": "TATANTATA",
            "X": "GATTACAGATTACAN"
        }
    }

The sequences are stored in a nested structure of Python dictionaries, each
key of the `sequence_data` dictionary represents the name or an otherwise unique
identifier for a particular sample (e.g. "my_sample", "my_other_sample"), the
value is a dictionary whose own keys represent chromosome names or numbers [#]_
and the corresponding values are the sequences themselves as a string [#]_.
Regardless of how the chromosomes are identified, they must match across samples
if one wishes to make comparisons across samples.

.. [#] Goldilocks has no preference for use of numbers or strings for chromosome names but it would be sensible to use numbers where possible for cases where you might wish to sort by chromosome.
.. [#] In future it is planned that sequences may be further nested in a dictionary to fully support polyploid species.


Conducting a Census
-------------------

Once you have organised your sequence data in to the appropriate structure, you
may conduct the census with Goldilocks by passing your strategy (e.g. NucleotideCounterStrategy)
and sequence data to the imported :class:`goldilocks.goldilocks.Goldilocks` class: ::

    g = Goldilocks(NucleotideCounterStrategy(["N"]), sequence_data, length=3, stride=1)

Make sure you add the brackets after the name of the imported strategy, this
'creates' a usuable strategy for Goldilocks to work with.

For each chromosome (i.e. 'one', 'X' and 2) Goldilocks will split each sequence
from the corresponding chromosome in each of the two example samples in to triplets
of bases (as our specified region length is 3) with an offset of 1 (as our stride is 1).
For example, chromosome `"one"` of `"my_sample"` will be split as follows: ::

    CAT
     ATC
      TCA
       CAN
        ANC
         NCA
          CAT

In our example, the NucleotideCounterStrategy will then count the number of N bases that
appear in each split, for each sample, for each chromosome.


Getting the Regions
-------------------

Once the census is complete, you can extract all of the censused regions directly
from your Goldilocks object. The example below demonstrates the format of the
returned regions dictionary for the example data above: ::

    > g.regions
    {
        0: {
            'chr': 2,
            'ichr': 0,
            'pos_end': 3,
            'pos_start': 1,
            'group_counts': {
                'my_sample': {'default': 2},
                'my_other_sample': {'default': 1},
                'total': {'default': 3}
            },
        }

        ...

        27: {
            'chr': 'one',
            'ichr': 6,
            'pos_end': 9,
            'pos_start': 7,
            'group_counts': {
                'my_sample': {'default': 0},
                'my_other_sample': {'default': 0},
                'total': {'default': 0}
            },
        }
    }


The returned structure is a dictionary whose keys represent the `id` of each region,
with values corresponding to a dictionary of metadata for that particular `id`.
The `id` is assigned incrementally (starting at 0) as each region is encountered
by Goldilocks during the census and isn't particularly important.

Each region dictionary has the following metadata [#]_:

============    =====
Key             Value
============    =====
id              A unique id assigned to the region by Goldilocks
chr             The chromosome the region appeared on (as found in the input data)
ichr            This region is the `ichr-th` to appear on this chromosome (0-indexed)
pos_start       The 1-indexed base of the sequence where the region begins (inclusive)
pos_end         The 1-indexed base of the sequence where the region ends (inclusive)
============    =====

.. [#] Goldilocks used to feature a group_counts dictionary as part of the region
       metadata as shown in the example above, this was removed as it duplicated
       data stored in the group_counts variable in the Goldilocks object needlessly.
       It has not been removed in the example output above as it helps explain
       what regions represent.


In the example output above, the first (0th) censused region appears on
chromosome 2 [#]_ and includes bases 1-3. It is the first (0th) region to appear on this
chromosome and over those three bases, the corresponding subsequence for `"my_sample"`
contained 2 N bases and the corresponding subsequence for `"my_other_sample"` contained
1. In total, over both samples, on chromosome 2, over bases 1-3, 3 N bases appeared.

The last region, region 27 (28th) appears on chromosome `"one"` [#]_ and includes
bases 7-9. It is the seventh (6th by 0-index) found on this chromosome and over
those three bases neither of the two samples contained an N base.

.. [#] As numbers are ordered before strings like "one" and "X" in Python.
.. [#] As "X" is ordered before "one" in Python.


Sorting Regions
---------------

Following a census, Goldilocks allows you to sort the regions found by four
mathematical operations: `max`, `min`, `mean` and `median`. ::

    regions_max = g.query("max")
    regions_min = g.query("min")
    regions_mean = g.query("mean")
    regions_median = g.query("median")

The data is returned in a special list:, a :class:`goldilocks.goldilocks.CandidateList`
which defines a table-based representation should a user wish to print the list: ::

    > print(regions_max)
    ID    VAL     CHR     POSITIONS (INC.)
    0       {'default': 3}  2                1 -          3
    2       {'default': 3}  2                3 -          5
    4       {'default': 3}  2                5 -          7
    6       {'default': 3}  2                7 -          9
    1       {'default': 2}  2                2 -          4
    ...
    18      {'default': 0}  X               11 -         13
    19      {'default': 0}  X               12 -         14
    21      {'default': 0}  one              1 -          3
    22      {'default': 0}  one              2 -          4
    27      {'default': 0}  one              7 -          9


Note the regions in the `regions_max` CandidateList are now sorted by the number
of N bases that appeared. Ties are currently resolved by the region that was seen
first (has the lowest `id`).


Full Example
------------

Census an example sequence for appearance of 'N' bases: ::

    from goldilocks.goldilocks import Goldilocks
    from goldilocks.strategies import NucleotideCounterStrategy

    sequence_data = {
        "my_sample": {
            2: "NANANANANA",
            "one": "CATCANCAT",
            "X": "GATTACAGATTACAN"
        },
        "my_other_sample": {
            2: "GANGANGAN",
            "one": "TATANTATA",
            "X": "GATTACAGATTACAN"
        }
    }

    g = Goldilocks(NucleotideCounterStrategy(["N"]), sequence_data, length=3, stride=1)

    regions_max_n_bases = g.query("max")
    regions_min_n_bases = g.query("min")
    regions_median_n_bases = g.query("min")
    regions_mean_n_bases = g.query("min")

