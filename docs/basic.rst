===================
Basic Package Usage
===================

The following example assumes basic Python programming experience (and
that you have installed Goldilocks), skip to the
end if you think you know what you're doing.

Importing
---------

To use Goldilocks you will need to import the :class:`goldilocks.goldilocks.Goldilocks`
class and your desired census strategy (e.g. NucleotideCounterStrategy) from
:mod:`goldilocks.strategies` to your script: ::

    from goldilocks import Goldilocks
    from goldilocks.strategies import NucleotideCounterStrategy


Providing Sequence Data as Dictionary
-------------------------------------

If you do not have FASTA files, the :class:`goldilocks.goldilocks.Goldilocks` class
allows you to provide sequence data in the following structure: ::

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


Providing Sequence Data as FASTA
--------------------------------
If your sequences are in FASTA format, you must first index them with `samtools faidx`,
then for each sample, pass the path to the index to Goldilocks in the following
structure: ::

    sequence_data = {
            "my_sequence": {"file": "/path/to/fastaidx/1.fa.fai"},
            "my_other_sequence": {"file": "/path/to/fastaidx/2.fa.fai"},
            "my_other_other_sequence": {"file": "/path/to/fastaidx/3.fa.fai"},
    }

When supplying sequences in this format, note the following:

    * `is_faidx=True` must be passed to the Goldilocks constructor (see below),
    * It is assumed that the FASTA will be in the same directory with the same name as its index, just without the ".fai" extension,
    * The key in the sequence data dictionary for each sample, must be `file`,
    * The `i`-th sequence in each FASTA will be censused together, thus the order in which your sequences appear matters.

It is anticipated in future these assumptions will be circumvented by additional
options to the Goldilocks constructor.

To specify the `is_faidx` argument, call the constructor like so: ::

    g = Goldilocks(NucleotideCounterStrategy(["N"]), sequence_data, length=3, stride=1, is_faidx=True)

Now Goldilocks will know to expect to open these `file` values as FASTA indexes,
not sequence data!

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

    g_max = g.query("max")
    g_min = g.query("min")
    g_mean = g.query("mean")
    g_median = g.query("median")

The result of a query is the original :class:`goldilocks.goldilocks.Goldilocks` object
with masked and sorted internal data. You can view a table-based representation
of the regions with :func:`goldilocks.goldilocks.Goldilocks.export_meta`. ::

    > g_max.export_meta(sep='\t', group="total")
    [NOTE] Filtering values between 0.00 and 3.00 (inclusive)
    [NOTE] 28 processed, 28 match search criteria, 0 excluded, 0 limit
    chr     pos_start       pos_end total_default
    2       1       3       3.0
    2       3       5       3.0
    2       5       7       3.0
    2       7       9       3.0
    2       2       4       2.0
    2       4       6       2.0
    2       6       8       2.0
    2       8       10      2.0
    X       13      15      2.0
    one     4       6       2.0
    one     5       7       2.0
    one     3       5       1.0
    one     6       8       1.0
    X       1       3       0.0
    X       2       4       0.0
    X       3       5       0.0
    X       4       6       0.0
    X       5       7       0.0
    X       6       8       0.0
    X       7       9       0.0
    X       8       10      0.0
    X       9       11      0.0
    X       10      12      0.0
    X       11      13      0.0
    X       12      14      0.0
    one     1       3       0.0
    one     2       4       0.0
    one     7       9       0.0

Note the regions in `g_max` are now sorted by the number
of N bases that appeared. Ties are currently resolved by the region that was seen
first (has the lowest `id`).

Setting Number of Processes
---------------------------

Goldilocks supports multiprocessing and can spawn some number of additional processes to perform
the census steps before aggregating all the region counters and answering queries.
To specify the number of processes Goldilocks should use, specify a `processes`
argument to the constructor: ::

    g = Goldilocks(NucleotideCounterStrategy(["N"]), sequence_data, length=3, stride=1, processes=4)

Full Example
------------

Census an example sequence for appearance of 'N' bases: ::

    from goldilocks import Goldilocks
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

    g = Goldilocks(NucleotideCounterStrategy(["N"]), sequence_data, length=3, stride=1, processes=4)

    g_max_n_bases = g.query("max")
    g_min_n_bases = g.query("min")
    g_median_n_bases = g.query("median")
    g_mean_n_bases = g.query("mean")

