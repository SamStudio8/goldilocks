==================
Command Line Usage
==================

Goldilocks is also packaged with a basic command line tool to demonstrate
some of its capabilities and to provide access to base functionality without
requiring users to author a script of their own. For more complicated queries,
you'll need to import Goldilocks as a package to a script of your own.
But for simple use-cases the tool might be enough for you.

Usage
#####

Goldilocks is invoked as follows: ::

    goldilocks <strategy> <sort-op> [--tracks TRACK1 [TRACK2 ...]] -l LENGTH -s STRIDE [-@ THREADS] FAIDX1 [FAIDX2 ...]

Where a strategy is a census strategy listed as available... ::

    $ goldilocks list
    Available Strategies
      * gc
      * ref
      * motif
      * nuc

...and a sort operation is one of:

* `max`
* `min`
* `mean`
* `median`
* `none`


Example
#######

Tabulate all regions and their associated counts of nucleotides A, C, G, T and N.
Window size 100Kbp, overlap 50Kbp. Census will spawn 4 processes. Regions in
table will be sorted by co-ordinate: ::

    goldilocks nuc none --tracks A C G T N -l 100000 -s 50000 -@ 4 /store/ref/hs37d5.fa.fai

Tabulate all regions and their associated GC-content. Same parameters as previous
example but table will be sorted by maximum GC-content descending: ::

    goldilocks gc max -l 100000 -s 50000 -@ 4 /store/ref/hs37d5.fa.fai

