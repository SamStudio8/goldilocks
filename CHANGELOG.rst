History
=======

0.0.80 (2015-08-10)
-------------------
* Added multiprocessing capabilities during census step.
* Added a simple command line interface.
* Removed prepare-evaluate paradigm from strategies and now perform counts
  directly on input data in one step.
* Skip slides (and set all counts to 0) if their `end_pos` falls outside of
  the region on that particular genome's chromosome/contig.
* Rename `KMerCounterStrategy` to `MotifCounterStrategy`
* Fixed bug causing `use_and` to not work as expected for chromosomes not
  explicitly listed in the `exceptions` dict when also using `use_chrom`.
* Support use of FASTA files which must be supplied with a `samtools faidx` style index.
* Stopped supporting Python 3 due to incompatability with `buffer` and `memoryview`.
* Prevent `query` from deep copying itself on return. Note this means that a query
  will alter the original Goldilocks object.
* Now using a 3D numpy matrix to store counters with memory shared to
  support multiprocessing during census.
* Removed `StrategyValue` as these cannot be stored in shared memory. This makes
  ratio-based strategies a bit of a hack currently (but still work...)
* tldr; Goldilocks is at least 2-4x faster than previously, even without multiprocessing

0.0.71 (2015-07-11)
-------------------
* Officially add MIT license to repository.
* Deprecate `_filter`.
* Update and tidy `examples.py`.
* `is_seq` argument to initialisation removed and replaced with `is_pos`.
* Use `is_pos` to indicate the expected input is positional, not sequence.
* Force use of `PositionCounterStrategy` when `is_pos` is True.
* Sequence data now read in to 0-indexed arrays to avoid the overhead of string
    re-allocation by having to append a padding character to the beginning of very
    long strings.
* Region metadata continues to use 1-indexed positions for user output.
* `VariantCounterStrategy` now `PositionCounterStrategy`.
* `PositionCounterStrategy` expects 1-indexed lists of positions;
    `prepare` populates the listed locations with 1 and then `evaluate`
    returns the sum as before.
* `test_regression2` updated to account for converting 1-index to 0-index when
    manually handling the sequence for expected results.
* `query` accepts `gmax` and `gmin` arguments to filter candidate regions by
  the group-track value.
* `CandidateList` removed and replaced with simply returning a new `Goldilocks`.

0.0.6 (2015-06-23)
------------------
* `Goldilocks.sorted_regions` stores a list of region ids to represent the result
  of a sorting operation following a call to `query`.
* Regions in `Goldilocks.regions` now always have a copy of their "id" as a key.
* `__check_exclusions` now accepts a `group` and `track` for more complex
  exclusion-based operations.
* `region_group_lte` and `region_group_gte` added to usable exclusion fields to
  remove regions where the value of the desired group/track combination is
  less/greater than or equal to the value of the group/track set by the
  current `query`.
* `query` now returns a new `Goldilocks` instance, rather than a `CandidateList`.
* `Goldilocks.candidates` property now allows access to regions, this property
  will maintain the order of `sorted_regions` if it has one.
* `export_meta` now allows `group=None`
* `CandidateList` class deleted.
* Test data that is no longer used has been deleted.
* Scripts for generating test data added to `test_gen/` directory.
* Tests updated to reflect the fact `CandidateList` lists are no longer returned
  by `query`.
* `_filter` is to be deprecated in favour of `query` by 0.0.7

Beta (2014-10-08)
---------------------
* Massively updated! Compatability with previous versions very broken.
* Software retrofitted to be much more flexible to support a wider range of problems.

0.0.2 (2014-08-18)
---------------------

* Remove incompatible use of `print`

0.0.1 (2014-08-18)
---------------------

* Initial package
