=========
Exporting
=========

Goldilocks provides functions for the exporting of all censused regions metadata
or for filtered regions resulting from a query. The examples below follow
on from the basic usage instructions earlier in the documentation.

Census Data
-----------

For a given sample one may export basic metadata for all regions that included
sequence data from that particular sample. The header is as follows:

====================    =====
Key                     Value
====================    =====
id                      A unique id assigned to the region by Goldilocks

track1                  The value for the region as calculated by the strategy used.
                        By default if a list of tracks is not specified when the
                        strategy is created, there will be just one track named
                        'default'.
                        For the majority of 'basic' strategies this will be the case.

[track2 ... trackN]     Optional further fields will appear for additional tracks,
                        the column header will feature the name of the track.
                        For example, a k-mer counting strategy would feature a
                        column for each k-mer specified to the strategy.

chr                     The chromosome the region appeared on (as found in the
                        input data)

pos_start               The 1-indexed base of the sequence where the
                        region begins (inclusive)

pos_end                 The 1-indexed base of the sequence where the region ends (inclusive)

====================    =====

Using the `my_sample` data: ::

    ...
    g.export_meta("my_sample", sep="\t")

    id      default chr     pos_start       pos_end
    0       2       2       1       3
    1       1       2       2       4
    2       2       2       3       5
    3       1       2       4       6
    4       2       2       5       7
    5       1       2       6       8
    6       2       2       7       9
    7       1       2       8       10
    8       0       X       1       3
    9       0       X       2       4
    10      0       X       3       5
    11      0       X       4       6
    12      0       X       5       7
    13      0       X       6       8
    14      0       X       7       9
    15      0       X       8       10
    16      0       X       9       11
    17      0       X       10      12
    18      0       X       11      13
    19      0       X       12      14
    20      1       X       13      15
    21      0       one     1       3
    22      0       one     2       4
    23      0       one     3       5
    24      1       one     4       6
    25      1       one     5       7
    26      1       one     6       8
    27      0       one     7       9

FASTA
-----

From any sorting or filtering operation on censused regions, a new Goldilocks
object is returned, providing function to output filtered sequence data to FASTA format.

Following on from the example introduced earlier, the example below shows the
subsequences of `my_sample` in the FASTA format, ordered by their appearance in
the filtered `candidates` list, from the highest number of 'N' bases, to the
lowest. ::

    ...
    candidates = g.query("max", group="my_sample")
    candidates.export_fasta("my_sample")

    >my_sample|Chr2|Pos1:3
    NAN
    >my_sample|Chr2|Pos3:5
    NAN
    >my_sample|Chr2|Pos5:7
    NAN
    >my_sample|Chr2|Pos7:9
    NAN
    >my_sample|Chr2|Pos2:4
    ANA
    >my_sample|Chr2|Pos4:6
    ANA
    >my_sample|Chr2|Pos6:8
    ANA
    >my_sample|Chr2|Pos8:10
    ANA
    >my_sample|ChrX|Pos13:15
    CAN
    >my_sample|Chrone|Pos4:6
    CAN
    >my_sample|Chrone|Pos5:7
    ANC
    >my_sample|Chrone|Pos6:8
    NCA
    >my_sample|ChrX|Pos1:3
    GAT
    >my_sample|ChrX|Pos2:4
    ATT
    >my_sample|ChrX|Pos3:5
    TTA
    >my_sample|ChrX|Pos4:6
    TAC
    >my_sample|ChrX|Pos5:7
    ACA
    >my_sample|ChrX|Pos6:8
    CAG
    >my_sample|ChrX|Pos7:9
    AGA
    >my_sample|ChrX|Pos8:10
    GAT
    >my_sample|ChrX|Pos9:11
    ATT
    >my_sample|ChrX|Pos10:12
    TTA
    >my_sample|ChrX|Pos11:13
    TAC
    >my_sample|ChrX|Pos12:14
    ACA
    >my_sample|Chrone|Pos1:3
    CAT
    >my_sample|Chrone|Pos2:4
    ATC
    >my_sample|Chrone|Pos3:5
    TCA
    >my_sample|Chrone|Pos7:9
    CAT
