History
=======

0.0.6 (2014-06-23)
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
