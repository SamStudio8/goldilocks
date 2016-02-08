========
Examples
========
The following includes some simple examples of what Goldilocks can be used for.

Example One
-----------

Read a pair of 1-indexed base position lists and output all regions falling
within 2 of the maximum count of positions in regions across both, in a table. ::

    from goldilocks import Goldilocks
    from goldilocks.strategies import PositionCounterStrategy
    data = {
            "my_positions": {
                1: [1,2,5,10,15,15,18,25,30,50,51,52,53,54,55,100]
            },
            "my_other_positions": {
                1: [1,3,5,7,9,12,15,21,25,51,53,59,91,92,93,95,99,100]
            }
    }
    g = Goldilocks(PositionCounterStrategy(), data, is_pos=True, length=10, stride=1)

    g.query("max", actual_distance=2).export_meta(sep="\t", group="total")

Example Two
-----------
Read a short sequence, census GC-ratio and output the top 5 regions as FASTA. ::

    from goldilocks import Goldilocks
    from goldilocks.strategies import GCRatioStrategy
    data = {
            "my_sequence": {
                1: "ACCGAGAGATTT"
            }
    }
    g = Goldilocks(GCRatioStrategy(), data, 3, 1)

    g.query("max", limit=5).export_fasta()

Example Three
-------------
Read a short sequence and census the appearance of the "AAA" and "CCC" motif.
Output a table of regions with the most occurrences of CCC (and at least one)
and another table of regions featuring the most appearances of both motifs.
Output only the maximum region (actual_distance = 0) displaying both motifs to
FASTA. ::

    from goldilocks import Goldilocks
    from goldilocks.strategies import MotifCounterStrategy
    data = {
            "my_sequence": {
                1: "CCCAAACCCGGGCCCGGGAGAAACCC"
            }
    }
    g = Goldilocks(MotifCounterStrategy(["AAA", "CCC"]), data, 9, 1)

    g.query("max", track="CCC", gmin=1).export_meta(sep="\t")
    g.query("max", group="total").export_meta(sep="\t", group="total", track="default")

    g.query("max", group="total", actual_distance=0).export_fasta()

Example Four
------------
Read two samples of three short chromosomes and search for 'N' nucleotides.
List and export a FASTA of regions that contain at least one N, sorted by number
of N's appearing across both samples. Below, an example of complex filtering. ::

    from goldilocks import Goldilocks
    from goldilocks.strategies import NucleotideCounterStrategy
    data = {
        "sample_one": {
            1: "ANAGGGANACAN",
            2: "ANAGGGANACAN",
            3: "ANANNNANACAN",
            4: "NNNNAANNAANN"
        },
        "sample_two": {
            1: "ANAGGGANACAN",
            2: "ANAGGGANACAN",
            3: "ANANNNANACAN",
            4: "NNNANNAANNAA"
        }
    }
    g = Goldilocks(NucleotideCounterStrategy(["N"]), data, 3, 1)

    g_max = g.query("max", gmin=1)
    g_max.export_meta(sep="\t")
    g_max.export_fasta()

    g.query("min",
            gmin = 1,
            exclusions={
                # Filter any region with a starting position <= 3 or >= 10
                "start_lte": 3,
                "start_gte": 10,

                # Filter any regions on Chr1
                1: {
                    "chr": True
                },

                # Filter NO regions on Chr2
                # NOTE: This also prevents the superexclusions above being applied.
                2: {
                    "chr": False
                },

                # Filter any region on Chr3 with an ending postion >= 9
                3: {
                    "start_lte": 5 # NOTE: This overrides the start_lte applied above
                }
            }, use_chrom=True).export_meta(sep="\t")

Example Five
------------
Read in four simple chromosomes from one sample and census the GC ratio.
Plot both a scatter plot of all censused regions over both of the provided
samples with position over the x-axis and value on the y-axis.
Produce a second plot drawing a panel with a line graph for each chromosome
with the same axes but data from one sample only.
For the combined result of both samples and chromosomes, organise the result
of the census for each region into desired bins and plot the result as a histogram.
Repeat the process for the my_sequence sample and produce a panelled histogram
for each chromosome. ::

    from goldilocks import Goldilocks
    from goldilocks.strategies import GCRatioStrategy
    data = {
        "my_sequence": {
            1: "ANAGGGANACANANAGGGANACANANAGGGANACANANAGGGANACANANAGGGACGCGCGCGGGGANACAN"*500,
            2: "ANAGGCGCGCNANAGGGANACGCGGGGCCCGACANANAGGGANACANANAGGGACGCGCGCGCGCCCGACAN"*500,
            3: "ANAGGCGCGCNANAGGGANACGCGGGGCCCGACANANAGGGANACANANAGGGACGCGCGCGCGCCCGACAN"*500,
            4: "GCGCGCGCGCGCGCGCGGGGGGGGGCGCCGCCNNNNNNNNNNNNNNNNGCGCGCGCGCGCGCGNNNNNNNNN"*500
        },
        "my_same_sequence": {
            1: "ANAGGGANACANANAGGGANACANANAGGGANACANANAGGGANACANANAGGGACGCGCGCGGGGANACAN"*500,
            2: "ANAGGCGCGCNANAGGGANACGCGGGGCCCGACANANAGGGANACANANAGGGACGCGCGCGCGCCCGACAN"*500,
            3: "ANAGGCGCGCNANAGGGANACGCGGGGCCCGACANANAGGGANACANANAGGGACGCGCGCGCGCCCGACAN"*500,
            4: "GCGCGCGCGCGCGCGCGGGGGGGGGCGCCGCCNNNNNNNNNNNNNNNNGCGCGCGCGCGCGCGNNNNNNNNN"*500
        }
    }
    g = Goldilocks(GCRatioStrategy(), data, 50, 10)

    g.plot()
    g.plot("my_sequence")
    g.profile(bins=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    g.profile("my_sequence", bins=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

Example Six
-----------
Read a set of simple chromosomes from two samples and tabulate the top 10% of
regions demonstrating the worst consensus to the given reference over both samples.
Plot the lack of consensus as line graphs for each chromosome, for each sample,
then over all chromosomes for all samples on one graph. ::

    from goldilocks import Goldilocks
    from goldilocks.strategies import ReferenceConsensusStrategy
    data = {
        "first_sample": {
            1: "NNNAANNNNNCCCCCNNNNNGGGGGNNNNNTTTTTNNNNNAAAAANNNNNCCCCCNNNNNGGGGGNNNNNTTTTTNNNNN",
            2: "NNNNNCCCCCNNNNNTTTTTNNNNNAAAAANNNNNGGGGGNNNNNCCCCCNNNNNTTTTTNNNNNAAAAANNNNNGGGGN"
        },
        "second_sample": {
            1: "NNNNNNNNNNCCCCCCCCCCNNNNNNNNNNTTTTTTTTTTNNNNNNNNNNCCCCCCCCCCNNNNNNNNNNTTTTTTTTTT",
            2: "NNCCCCCCCCNNNNNNNNNNAAAAAAAAAANNNNNNNNNNCCCCCCCCCCNNNNNNNNNNAAAAAAAAAANNNNNNNNNN"
        }
    }
    ref = {
        1: "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT",
        2: "CCCCCCCCCCTTTTTTTTTTAAAAAAAAAAGGGGGGGGGGCCCCCCCCCCTTTTTTTTTTAAAAAAAAAAGGGGGGGGGG"
    }

    g = Goldilocks(ReferenceConsensusStrategy(reference=ref, polarity=-1), data, stride=10, length=10)
    g.query("max", percentile_distance=10).export_meta(group="total", track="default")

    g.plot("first_sample")
    g.plot("second_sample")
    g.plot()

Example Seven
-------------
Read a pair of 1-indexed base position lists from two samples. Sort regions
with the least number of marked positions on Sample 1 and subsort by max marked
positions in Sample 2. ::

    from goldilocks import Goldilocks
    from goldilocks.strategies import PositionCounterStrategy
    data = {
            "my_positions": {
                1: [1,2,3,4,5,6,7,8,9,10,
                    11,13,15,17,19,
                    21,
                    31,39,
                    41]
            },
            "other_positions": {
                1: [21,22,23,24,25,26,27,28,
                    31,33,39,
                    41,42,43,44,45,46,47,48,49,50]
            }
    }
    g = Goldilocks(PositionCounterStrategy(), data, is_pos=True, length=10, stride=5)

    g.query("max", group="my_positions").query("max", group="other_positions").export_meta(sep="\t")
