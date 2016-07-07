==========
Goldilocks
==========

.. image:: https://badge.fury.io/py/goldilocks.png
    :target: http://badge.fury.io/py/goldilocks

.. image:: https://travis-ci.org/SamStudio8/goldilocks.png?branch=master
        :target: https://travis-ci.org/SamStudio8/goldilocks

.. image:: https://coveralls.io/repos/SamStudio8/goldilocks/badge.png?branch=master
        :target: https://coveralls.io/r/SamStudio8/goldilocks

.. image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Join the chat at https://gitter.im/SamStudio8/goldilocks
   :target: https://gitter.im/SamStudio8/goldilocks?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge

Locating genomic regions that are "just right".

* Documentation: http://goldilocks.readthedocs.org.


What is it?
-----------

**Goldilocks** is a Python package providing functionality for locating 'interesting'
genomic regions for some definition of 'interesting'. You can import it to your
scripts, pass it sequence data and search for subsequences that match some criteria
across one or more samples.

Goldilocks was developed to support our work in the investigation of quality
control for genetic sequencing. It was used to quickly locate
regions on the human genome that expressed a desired level of variability,
which were "just right" for later variant calling and comparison.

The package has since been made more flexible and can be used to find regions
of interest based on other criteria such as GC-content, density of target k-mers,
defined confidence metrics and missing nucleotides.


What can I use it for?
----------------------

Given some genetic sequences (from one or more samples, comprising of one or more
chromosomes), Goldilocks will shard each chromosome in to subsequences of a
desired size which may or may not overlap as required. For each chromosome from
each sample, each subsequence or 'region' is passed to the user's chosen strategy.

The strategy simply defines what is of interest to the user in a language that
Goldilocks can understand. Goldilocks is currently packaged with the following
strategies:

============================      ==================
Strategy                          Census Description
============================      ==================
GCRatioStrategy                   Calculate GC-ratio for subregions across the
                                  genome.
NucleotideCounterStrategy         Count given nucleotides for subregions across
                                  the genome.
MotifCounterStrategy              Search for one or more particular motifs of
                                  interest of any and varying size in subregions
                                  across the genome.
ReferenceConsensusStrategy        Calculate the (dis)similarity to a given
                                  reference across the genome.
PositionCounterStrategy           Given a list of base locations, calculate
                                  density of those locations over subregions
                                  across the genome.
============================      ==================

Once all regions have been 'censused', the results may be sorted by one of four
mathematical operations: `max`, `min`, `median` and `mean`. So you may be interested
in subregions of your sequence(s) that feature the most missing nucleotides, or
subregions that contain the mean or median number of SNPs or the lowest GC-ratio.


Why should I use it?
--------------------

Goldilocks is hardly the first tool capable of calculating GC-content across a
genome, or to find k-mers of interest, or SNP density, so why should you use it
as part of your bioinformatics pipeline?

Whilst not the first program to be able to conduct these tasks, it is the first
to be capable of doing them all together, sharing the same interfaces. Every strategy
can quickly be swapped with another by changing one line of your code. Every strategy
returns regions in the same format and so you need not waste time munging data to
fit the rest of your pipeline.

Strategies are also customisable and extendable, those even vaguely familiar with
Python should be able to construct a strategy to meet their requirements.

Goldilocks is maintained, documented and tested, rather than that hacky perl
script that you inherited years ago from somebody who has now left your lab.


Requirements
------------
To use;

* numpy
* matplotlib (for plotting)

To test;

* tox
* pytest

For coverage;

* nose
* python-coveralls

Installation
------------

::

    $ pip install goldilocks


Citation
--------

Please cite us so we can continue to make useful software! ::

    Nicholls, S. M., Clare, A., & Randall, J. C. (2016). Goldilocks: a tool for identifying genomic regions that are "just right." Bioinformatics (2016) 32 (13): 2047-2049. doi:10.1093/bioinformatics/btw116
    
::

    @article{Nicholls01072016,
        author = {Nicholls, Samuel M. and Clare, Amanda and Randall, Joshua C.}, 
        title = {Goldilocks: a tool for identifying genomic regions that are ‘just right’},
        volume = {32}, 
        number = {13}, 
        pages = {2047-2049}, 
        year = {2016}, 
        doi = {10.1093/bioinformatics/btw116}, 
        URL = {http://bioinformatics.oxfordjournals.org/content/32/13/2047.abstract}, 
        eprint = {http://bioinformatics.oxfordjournals.org/content/32/13/2047.full.pdf+html}, 
        journal = {Bioinformatics} 
    }

License
-------
Goldilocks is distributed under the MIT license, see LICENSE.
