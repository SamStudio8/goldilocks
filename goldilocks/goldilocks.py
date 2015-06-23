__author__ = "Sam Nicholls <sn8@sanger.ac.uk>"
__copyright__ = "Copyright (c) Sam Nicholls"
__version__ = "0.0.6"
__maintainer__ = "Sam Nicholls <sam@samnicholls.net>"

import numpy as np
from math import floor, ceil

import copy
import os

# TODO Generate database of regions with stats... SQL/SQLite
#      - Probably more of a wrapper script than core-functionality: goldib
# TODO Support more interesting sequence formats? FASTQ reading?
#      - Read FASTQ quality data but still output sequences (read both Q and SEQ)
# TODO Replace 'group' nonclementure with 'sample'?
class Goldilocks(object):
    """Facade class responsible for conducting a census of genomic regions.

    Given sequence data and a :mod:`goldilocks.strategies` search strategy,
    Goldilocks is able to census regions along multiple genomes of a desired
    length and overlap and provides an interface to query results for a given
    criteria.

    .. note:: **_filter deprecated in 0.0.53**

        `_filter` will be removed in Goldilocks 1.0.0, it is replaced by the
        more suitably named `query`. Underscores are traditionally used for
        private class members and `_filter` was merely named to prevent
        confusion with the built-in `filter`.
        Until 1.0, `_filter` will pass all of its arguments to `query`.

    Parameters
    ----------
    strategy : Strategy object
        An instantiated :mod:`goldilocks.strategies` search strategy.

    data : dict{str, dict{[str|int], str}}
        Data on which to conduct a census, in a nested dict.

        Top level keys represent individual samples whose value is another
        dict using str or int chromosome keys and sequence data as values.
        As an example: ::

            "my_sample": {
                "chrom_name_or_number": "SEQUENCE",
            },
            "my_other_sample": {
                "chrom_name_or_number": "SEQUENCE",
            }

    length : int
        Desired region length, all censused regions will be of this many bases.
        If a region will end such that it would end beyond the end of a chromosome,
        it will not be added to the census.

    stride : int
        Number of bases to add to the start of the last region before the
        start of the next. If LENGTH==STRIDE, there will be no overlap and
        regions will begin on the base following where the previous region ended.

    is_seq : boolean, optional(default=True)
        Whether or not the data stored in `data` is sequence data.
        If `is_seq` is False, Goldilocks will expect a list of base positions.

    Attributes
    ----------
    strategy : Strategy object
        The desired :mod:`goldilocks.strategies` search strategy as selected by
        the user upon initialisation. Goldilocks will use the exposed `prepare`
        and `evaluate` functions of `stategy` to conduct the census.

    data : dict{str, dict{[str|int], str}}
        Data on which to conduct a census, in a nested dict.

    length : int
        Desired region length, as provided by the user.

    stride : int
        Desired region stride, as provided by the user.

    chr_max_len : dict{str, int}
        Maps names of chromosomes to the largest size encountered for that
        chromosome across all samples

    group_counts : dict{str, dict{str, list{int}}}
        Each group contains a dictionary of track-counter lists.
        For each group-track pair, a list stores the values returned from
        strategy evaluation for each subregion encountered by the census
        function. Each value is appended to the relevant group-track list.

        When a strategy uses multiple groups and tracks, the region id is
        used to update the corresponding elements in these lists.

        Once the census is complete the data stored in these counters are
        used for calculating a target value such as the maximum, minimum,
        mean or median.

    group_buckets : dict{str, dict{str, dict{[int|float], list{int}}}}
        Each group contains a dictionary of track-bucket dicts.
        For each group-track, a dict maps values returned from strategy evaluation
        to a list of region ids that was evaluated to that value.

        In a very basic example where a census is conducted for 'A' nucleotides
        over one sample (group) which features one chromosome of length 16: ::

            1|AAAA..AA.AA.AAAA|16

        With a length of 4 and a stride of 4 (ie. an overlap of 0):

            ===   =====   ========  ===   =====
            ID    Start   Sequence  End   Value
            ===   =====   ========  ===   =====
            0     1       AAAA      4     4
            1     5       ..AA      8     2
            2     9       .AA.      12    2
            3     13      AAAA      16    4
            ===   =====   ========  ===   =====

        The buckets would be organised as thus: ::

            \ 2 /\ 4 /
              |    |
              |    > [0,3]
              |
              > [1,2]

        Once the desired 'target' value has been calculated (max, min, mean
        or median), these buckets are used to selected regions (by their ID)
        that fall inside the desired distance from the target without
        requiring iteration over all censused regions again.

    regions : dict{int, dict{str, int}}
        A dict mapping automatically assigned ascending (from 0) integer ids
        to censused region metadata including the following keys:

            =========    ===============================================
            Key          Value
            =========    ===============================================
            chr          Chromosome on which this region appears
            ichr         Region is the i'th to appear on this chromosome
            pos_start    1-indexed base this region starts on (incl.)
            pos_end      1-indexed base this region ends on (incl.)
            =========    ===============================================

        The region id can be used to access the corresponding counter
        information from the lists stored in group_counts.

        These ids are also the same ids saved in relevant group_buckets.

    sorted_regions : list{int}
        A list of region ids representing the result of a sorting operation
        after a call to `query`.

    MULTI_TRACKED : boolean
        Whether or not the selected search strategy is using more than
        one 'default' track. If this is true, the group_counts and group_buckets
        attributes will hold extra keys to summarise data over the various
        samples and tracks.

    Raises
    ------
    ValueError
        If either `length` or `stride` are less than one.

    Notes
    -----
    * Setting `is_seq` to False may currently be incompatible with some strategy types.

    """
    def __init__(self, strategy, data, length, stride, is_seq=True):

        self.strategy = strategy

        self.chr_max_len = {}

        self.groups = {}

        self.group_counts = {"total": {}}
        self.group_buckets = {"total": {}}

        self.regions = {}
        self.sorted_regions = []
        self.target = None

        # Read data
        self.groups = data
        self.max_chr_max_len = None

        # Is this strategy multi-tracked?
        self.MULTI_TRACKED = False
        if (len(self.strategy.TRACKS) > 1
                or (len(self.strategy.TRACKS) == 1 and "default" not in self.strategy.TRACKS)):
            #TODO Warn if trying to use a track named strategy
            self.MULTI_TRACKED = True

        for group in self.groups:
            #TODO Catch duplicates etc..
            for chrom in self.groups[group]:
                # Remember to exclude 0-index
                if is_seq:
                    # TODO This is awful.
                    # Prepend _ to sequence to force 1-index
                    if self.groups[group][chrom][0] != "_":
                        self.groups[group][chrom] = "_" + self.groups[group][chrom]

                    len_current_seq = len(self.groups[group][chrom]) - 1
                else:
                    len_current_seq = sorted(self.groups[group][chrom])[-1]

                #TODO Should we not be looking for min length?
                if chrom not in self.chr_max_len:
                    self.chr_max_len[chrom] = len_current_seq

                if len_current_seq > self.chr_max_len[chrom]:
                    self.chr_max_len[chrom] = len_current_seq

        # Initialise group-track counts and buckets
        for group in self.groups:
            self.group_buckets[group] = {}
            self.group_counts[group] = {}

            # Initialise storage for tracks
            for track in self.strategy.TRACKS:
                self.group_buckets[group][track] = {}

                self.group_counts[group][track] = []
                self.group_counts["total"][track] = []

            # Populate additional group counters if using more than just the 'default' track
            if self.MULTI_TRACKED:
                self.group_counts[group]["default"] = []

        if self.MULTI_TRACKED:
            self.group_counts["total"]["default"] = []

        # Ensure stride and length are valid (>1)
        if stride < 1:
            raise ValueError("[FAIL] Stride must be at least 1 base wide.")
        self.STRIDE = stride

        if length < 1:
            raise ValueError("[FAIL] Length must be at least 1 base wide.")
        self.LENGTH = length

        # Automatically conduct census
        self.census()

    def __bucketize(self, scores_to_bucket):
        buckets = {}
        for region_id, total in enumerate(scores_to_bucket):
            if total not in buckets:
                buckets[total] = []
            buckets[total].append(region_id)
        return buckets

    def __load_chromosome(self, arr, data, track):
        return self.strategy.prepare(arr, data, track)

    def census(self):
        """Conduct a census of genomic subregions of a given size and overlap over
        chromosomes identified in each submitted sample. For each chromosome, each
        sample is loaded and split in to regions of the correct size and overlap
        which are then processed and evaluated by the desired strategy."""

        # We'd like to try and do things in order, it can be useful to a user
        # who'd like to extract and plot the lists stored in group_counts directly.
        # Python3 dislikes us sorting mixed lists with ints and strings, so they're
        # ordered seperately and concatenated together.
        #NOTE Using 'type' is considered a little naughty but I'm not expecting
        #     use of subclassed ints, longs or anything overly complex here, but
        #     just in case, the code falls back to the chr_max_len item list.
        #     (As these objects would fall in to chr_str and a TypeError would
        #     be encountered due to the unorderable types)
        try:
            chr_num = [chrom for chrom in self.chr_max_len.items() if type(chrom[0])==int]
            chr_str = [chrom for chrom in self.chr_max_len.items() if type(chrom[0])!=int]

            chroms = sorted(chr_num)
            chroms.extend(sorted(chr_str))
        except TypeError:
            chroms = self.chr_max_len.items()

        regions = {}
        region_i = 0
        for chrno, size in chroms:
            chros = {}
            for group in self.groups:
                chros[group] = {}
                for track in self.strategy.TRACKS:
                    chro = np.zeros(size+1, np.int8)
                    chros[group][track] = self.__load_chromosome(chro, self.groups[group][chrno], track)

            print("[SRCH] Chr:%s" % str(chrno))
            # Ignore 0 position
            for i, region_s in enumerate(range(1, size+1-self.LENGTH+1, self.STRIDE)):
                region_e = region_s + self.LENGTH - 1
                regions[region_i] = {
                    "id": region_i,
                    "ichr": i,
                    "chr": chrno,
                    "pos_start": region_s,
                    "pos_end": region_e
                }

                for group in self.groups:
                    for track in self.strategy.TRACKS:
                        # Evaluate the prepared region using the selected strategy and track
                        value = self.strategy.evaluate(chros[group][track][region_s:region_e+1], track=track)

                        # TODO Should we be ignoring these regions if they are empty?
                        # TODO Config option to include 0 in flter metrics
                        if value not in self.group_buckets[group][track]:
                            # Add this particular number of variants as a bucket
                            self.group_buckets[group][track][value] = []

                        # Add the region id to the bucket
                        self.group_buckets[group][track][value].append(region_i)

                        self.group_counts[group][track].append(value)

                        try:
                            self.group_counts["total"][track][region_i] += value
                        except IndexError:
                            self.group_counts["total"][track].append(value)

                        if self.MULTI_TRACKED:
                            try:
                                self.group_counts["total"]["default"][region_i] += value
                            except IndexError:
                                self.group_counts["total"]["default"].append(value)

                            try:
                                self.group_counts[group]["default"][region_i] += value
                            except IndexError:
                                self.group_counts[group]["default"].append(value)

                region_i += 1

        # Populate total group_buckets now census is complete
        self.group_buckets["total"] = {}
        super_totals = []

        ggroups = list(self.groups.keys())
        ggroups.append("total")

        for count_group in ggroups:
            group_totals = []
            for track in self.strategy.TRACKS:
                totals = []
                for region_id, value in enumerate(self.group_counts[count_group][track]):
                    try:
                        totals[region_id] += value
                    except IndexError:
                        totals.append(value)

                    if count_group == "total":
                    # Total track already stores sums, don't want to count
                    # everything twice by using tracks in the non-total group
                        try:
                            super_totals[region_id] += value
                        except IndexError:
                            super_totals.append(value)
                    else:
                        try:
                            group_totals[region_id] += value
                        except IndexError:
                            group_totals.append(value)

                self.group_buckets[count_group][track] = self.__bucketize(totals)

            self.group_buckets[count_group]["default"] = self.__bucketize(group_totals)

        # Populate super total-default group-track which sums the totals across
        # all tracks in all groups
        self.group_buckets["total"]["default"] = self.__bucketize(super_totals)

        self.regions = regions

    #TODO Check that doubling the window size for max and min works as expected
    #TODO Is changing the behaviour (re: window) for different math ops a bad idea...
    def __apply_filter_func(self, func, lower_window, upper_window, group, actual, track):

        target = None
        if func.lower() == "median":
            if actual:
                q_low  = np.percentile(np.asarray(self.group_counts[group][track]), 50) - lower_window
                q_high = np.percentile(np.asarray(self.group_counts[group][track]), 50) + upper_window
            else:
                q_low  = np.percentile(np.asarray(self.group_counts[group][track]), 50 - lower_window)
                q_high = np.percentile(np.asarray(self.group_counts[group][track]), 50 + upper_window)
            target = np.percentile(np.asarray(self.group_counts[group][track]), 50)
        elif func.lower() == "mean":
            mean = np.mean(np.asarray(self.group_counts[group][track]))
            if actual:
                q_low  = mean - lower_window
                q_high = mean + upper_window
            else:
                # A crude (but probably 'close enough') calculation for the
                # mean's percentile standing.
                track_scores = np.asarray(self.group_counts[group][track])

                # Center scores around the mean
                track_scores = track_scores - mean

                # Divide number of scores <= 0 (the mean) by n
                #TODO Should probably be interpolating for cases where the mean
                #     is not contained in track_scores...?
                mean_percentile = (len(track_scores[track_scores <= 0]) / float(len(track_scores))) * 100

                q_low  = np.percentile(np.asarray(self.group_counts[group][track]), mean_percentile - lower_window)
                q_high = np.percentile(np.asarray(self.group_counts[group][track]), mean_percentile + upper_window)

            target = mean
        elif func.lower() == "max":
            q_high = np.percentile(np.asarray(self.group_counts[group][track]), 100)
            if actual:
                q_low  = np.percentile(np.asarray(self.group_counts[group][track]), 100) - lower_window
            else:
                q_low  = np.percentile(np.asarray(self.group_counts[group][track]), 100 - lower_window)
            target = q_high
        elif func.lower() == "min":
            q_low = np.percentile(np.asarray(self.group_counts[group][track]), 0)
            if actual:
                q_high = np.percentile(np.asarray(self.group_counts[group][track]), 0) + upper_window
            else:
                q_high = np.percentile(np.asarray(self.group_counts[group][track]), 0 + upper_window)
            target = q_low
        else:
            raise TypeError("[FAIL] Invalid sorting function: '%s'" % func.lower())
        return float(q_low), float(q_high), float(target)


    def __check_exclusions(self, exclusions, region_dict, group, track, use_and=False, use_chrom=False):
        if exclusions is None or len(exclusions) == 0:
            return False

        def __exclude_chro(region_dict, chr_list_or_bool):
            try:
                if region_dict["chr"] in chr_list_or_bool:
                    return True
            except TypeError:
                # It's probably a bool not a list (it could be an int if users
                # have failed to read the documentation...)
                if chr_list_or_bool:
                    return True
            return False

        def __exclude_start(region_dict, operand, position):
            if operand < 0:
                return region_dict["pos_start"] <= position
            elif operand > 0:
                return region_dict["pos_start"] >= position
            return False

        def __exclude_end(region_dict, operand, position):
            if operand < 0:
                return region_dict["pos_end"] <= position
            elif operand > 0:
                return region_dict["pos_end"] >= position
            return False

        #TODO A bit pointless
        def __exclude_val(threshold, operand, other_val):
            if operand < 0:
                return other_val <= threshold
            elif operand > 0:
                return other_val >= threshold
            return False

        # Chrom specific exclusions override overall (might be unexpected)
        if use_chrom:
            if region_dict["chr"] in exclusions:
                to_apply = exclusions[region_dict["chr"]]
            else:
                # No exclusions to apply to this region
                return False
        else:
            if region_dict["chr"] in exclusions:
                print("[WARN] Exclusions dictionary appears to contain the name of a chromosome. Did you forget to set use_chrom=True?")
            to_apply = exclusions

        for name in to_apply:
            #TODO Could probably improve with a dict of funcs...
            ret = False
            if name == "chr":
                ret = __exclude_chro(region_dict, to_apply["chr"])
            elif name == "start_lte":
                ret = __exclude_start(region_dict, -1, to_apply["start_lte"])
            elif name == "start_gte":
                ret = __exclude_start(region_dict, 1, to_apply["start_gte"])
            elif name == "end_lte":
                ret = __exclude_end(region_dict, -1, to_apply["end_lte"])
            elif name == "end_gte":
                ret = __exclude_end(region_dict, 1, to_apply["end_gte"])
            elif name == "region_group_lte":
                ret = __exclude_val(self.group_counts[group][track][region_dict["id"]], -1, self.group_counts[to_apply["region_group_lte"]][track][region_dict["id"]])
            elif name == "region_group_gte":
                ret = __exclude_val(self.group_counts[group][track][region_dict["id"]], 1, self.group_counts[to_apply["region_group_lte"]][track][region_dict["id"]])

            else:
                #TODO Better handling of invalid exclusion property
                print("[WARN] Attempted to exclude on invalid property '%s'" % name)

            if use_and:
                # Require all exclusions to be true...
                if not ret:
                    return False
            else:
                # If we're not waiting on all conditions, we can exclude on the first
                if ret:
                    return True

        if use_and:
            # If we didn't bail on a previous false, all conditions must be satisfied
            return True
        return False

    def query(self, func="median", track="default", actual_distance=None, percentile_distance=None,
            direction=0, group="total", limit=0, exclusions=None, use_and=False, use_chrom=False):
        """Query the Goldilocks census to retrieve regions that meet given criteria.

        Parameters
        ----------
        func : str, options={"max", "min", "mean", "median"}
            Sorting function to be applied when returning regions which meet
            the input criteria and from which to calculate distances for use
            with `actual_distance` or `percentile_distance`.

        group : str, optional(default="total")
            Sort and filter only using values evaluated by the strategy only from
            data in the given sample `group`. Whilst it is possible to census
            data from many different sample groups in the same census, you may
            only be interested in isolating regions on a particular sample.

            If a `group` is not provided, by default, the "total" `group` is used,
            which represents the aggregate (but not necessarily `sum`) of values
            seen at a region site across all samples.

            For example, the "total" `group` for a nucleotide counting strategy
            would contain the sum of all counted bases for all genomic
            sub-sequences that lie on the given region for each track in the
            strategy.

        track : str, optional(default="default")
            Sort and filter only using values evaluated by the strategy only from
            data in the given `track`. In a simple nucleotide counting example,
            whilst you can count for multiple bases in the same census, you may
            wish to query based on data from just 'N' bases.

            If a `track` is not provided, by default, the "default" `track` is used,
            which represents the aggregate (but not necessarily `sum`) of values
            seen at a region site across all tracks.

            For example, the "default" `track` for a nucleotide counting strategy
            would contain the sum of all counted bases over a given region for
            a given group.

            .. note::
                The "total" `group` contains a "default" `track` which for a simple
                nucleotide counting strategy, would hold the sum of all bases of
                interest seen across sub-sequences on all groups, over all tracks
                that lie on the given region.

            .. note::
                Ratio-based strategies that do not simply count instances of given
                bases or motifs etc. will be correctly weighted by use of
                :class:`goldilocks.strategies.StrategyValue` in aggregate groups
                and tracks. No special handling for these sort of strategies is required.

                However, if you are writing a custom strategy, be sure to return
                a :class:`goldilocks.strategies.StrategyValue` from your `evaluate`
                method.

        actual_distance : float, optional(default=None)
        percentile_distance : float, optional(default=None)
            Filter regions whose value as returned from the selected strategy
            falls outside the absolute or percentile distance from the target
            as calculated by `func`.

            `actual_distance` will filter regions whose difference from the
            target falls outside the given value.

            `percentile_distance` will filter regions whose value falls outside
            the given number of percentiles.

            When used with `direction` one may decide whether to look above,
            below or around the target.

            .. note::
                `actual_distance` and `percentile_distance` are mutually exclusive.

        direction : int, optional(default=0)
            When using `actual_distance` or `percentile_distance` one may select
            whether to select regions that appear within the desired distance
            above, below or around the target - as calculated by `func`.

            Any positive value will set the direction to "upper", any negative
            value will set the direction to "lower". By default `direction` is 0,
            which will search around the target.

            For example, to find the 25% of values that appear above the mean, set
            `percentile_distance` to 25, `func` to "mean" and `direction` to one.
            To find the 10% of values below the median, set `percentile_distance`
            to 10, `func` to "median" and `direction` to -1.

            To find regions within plus/minus 5.0 of the mean, set `func` to mean and
            `direction` to 0 and `actual_distance` to 10.

            .. note::
                If `func` is max or min, the direction will automatically
                be changed to +1 or -1, respectively - as it doesn't make sense to
                search "around" the maximum or minimum value.

        limit : int, optional(default=0)
            Maximum number of regions to return.
            By default, all regions that meet the specified criteria will be returned.

        exclusions : [dict{str, dict{str, [int|str|list]}} | dict{[int|str], dict{str, [int|str|list|boolean]}}], optional(default=None)
            A dict defining criteria with which to filter regions.

            The dict may be specified in two ways: keys can either be
            exclusion properties as found in the table below or match the names
            or numbers of chromosomes provided to the constructor of Goldilocks
            with dict values that specify exclusion properties to values.

            The former method will apply specified exclusions to all regions
            whereas the latter when `use_chrom` is set to True will apply
            exclusions to particular chromosomes.

            Currently the following excluding criteria are available:

            =========   ==============================================================
            Criterion   Purpose
            =========   ==============================================================
            start_lte   Region starts on 1-indexed base less than or equal to value
            start_gte   Region starts on 1-indexed base greater than or equal to value
            end_lte     Region ends on 1-indexed base less than or equal to value
            end_gte     Region ends on 1-indexed base greater than or equal to value
            chr         Region appears on chr in given list
            =========   ==============================================================

            Further information and examples on using these effectively can be
            found in the documentation on sorting and filtering.

        use_and : boolean, optional(default=False)
            A flag to indicate whether a region must meet all exclusion criteria
            defined in `exclusions` to be excluded. By default this is False and
            a region will be excluded if it meets one or more exclusion criteria.

        use_chrom : boolean, optional(default=False)
            A flag to indicate that the keys of the `exclusions` dict are chromosome
            identifiers and exclusion criteria within should be applied to particular
            chromosomes. By default it is assumed that keys in the `exclusions`
            dict are to be applied to all regions, regardless of the chromosome
            on which they appear.

            .. note::
                `use_chrom` can be used with `use_and`, all criteria in each block
                of chromosome specific exclusions must be met for a region on that
                chromosome to be excluded.

            .. note::
                Goldilocks will print a warning to stdout if it encounters the name of a
                chromosome in the `exclusions` dict without `use_chrom` being set
                to true, but will continue to complete the query anyway.

        Returns
        -------
        Goldilocks object : :class:`goldilocks.goldilocks.Goldilocks`
            Returns a new instance of Goldilocks where the regions are filtered
            and the result of any sorting operation is stored in sorted_regions.
            Sorts are descending from absolute distance to the target as calculated
            by `func`. The target is also stored as `target`.

        Raises
        ------
        TypeError
            When attempting to sort by an invalid `func`.
        ValueError
            If attempting to filter both by `actual_distance` and `percentile_distance`.

        """
        #TODO Raise error when setting use_and or use_chrom without exclusions?
        #TODO Raise error when using direction without using a distance method?
        #     Probably more suitable to just raise a warning?

        distance = None
        actual = False

        if actual_distance is not None and percentile_distance is not None:
            raise ValueError("[FAIL] Cannot filter by both actual_distance and percentile_difference. Select one.")

        if actual_distance is not None:
            distance = actual_distance
            actual = True

        if percentile_distance is not None:
            distance = percentile_distance

        candidates = []
        if distance is not None:
            upper_window = float(distance)
            lower_window = float(distance)
            if func not in ["max", "min"]:
                upper_window = float(distance) / 2
                lower_window = float(distance) / 2

            #TODO ValueError when trying to use + direction on max and - direction on min?
            if direction > 0:
                upper_window = upper_window * 2
                lower_window = 0.0
            elif direction < 0:
                lower_window = lower_window * 2
                upper_window = 0.0

            q_low, q_high, target = self.__apply_filter_func(func, upper_window, lower_window, group, actual, track)
        else:
            #TODO Pretty ugly.
            q_low, q_high, target = self.__apply_filter_func(func, 0, 0, group, actual, track)
            q_low  = min(np.asarray(self.group_counts[group][track]))
            q_high = max(np.asarray(self.group_counts[group][track]))

        # For each "number of variants" bucket: which map the number of variants
        # seen in a region, to all regions that contained that number of variants
        print("[NOTE] Filtering values between %.2f and %.2f (inclusive)" % (floor(q_low), ceil(q_high)))

        for bucket in self.group_buckets[group][track]:
            if bucket >= floor(q_low) and bucket <= ceil(q_high):
                # Append all region data structures within the desired range
                # to the list of candidates for enrichment
                candidates += self.group_buckets[group][track][bucket]

        num_selected = 0
        num_excluded = 0
        num_total = 0

        to_delete = []
        sorted_regions = []
        for region in sorted(self.regions,
                    key=lambda x: (abs(self.group_counts[group][track][x] - target))):

            num_total += 1
            chosen = False
            if region in candidates:
                num_selected += 1
                if not self.__check_exclusions(exclusions, self.regions[region], group, track, use_and, use_chrom):
                    chosen = True
                    sorted_regions.append(region)
                else:
                    num_excluded += 1

            if not chosen:
                to_delete.append(region)

        # Return the top N elements if desired
        # TODO Report total, selected, selected-excluded and selected-filtered
        if limit > 0:
            for r in sorted_regions[limit:]:
                if r not in to_delete:
                    to_delete.append(r)
            sorted_regions = sorted_regions[0:limit]

        print("[NOTE] %d processed, %d match search criteria, %d excluded, %d limit" %
                (num_total, num_selected, num_excluded, limit))

        # TODO Pretty gross, it is probably worth brining back the CandidateList
        # object as a @property that provides the query function and then returning
        # frames from it for function chaining - rather than the Goldilocks instance...
        new_g = copy.deepcopy(self)
        for region in to_delete:
            del new_g.regions[region]

        print("[NOTE] %d regions pruned" % (len(to_delete)))
        new_g.sorted_regions = sorted_regions
        new_g.target = target
        return new_g

    #TODO Remove defaults
    def _filter(self, func="median", track="default", actual_distance=None, percentile_distance=None,
            direction=0, group="total", limit=0, exclusions=None, use_and=False, use_chrom=False): # pragma: no cover
        return self.query(func=func,
                          track=track,
                          actual_distance=actual_distance,
                          percentile_distance=percentile_distance,
                          direction=direction,
                          group=group,
                          limit=limit,
                          exclusions=exclusions,
                          use_and=use_and,
                          use_chrom=use_chrom)

    def plot(self, group=None, track="default", ylim=None, save_to=None, annotation=None): # pragma: no cover
        """Represent censused regions in a plot using matplotlib."""

        import matplotlib.pyplot as plt

        if group is None:
            group = "total"
            fig = plt.subplot(1,1,1)

            num_regions = len(self.regions)
            num_counts = [self.group_counts[group][track][x] for x in sorted(self.regions)]

            max_val = max(num_counts)
            plt.scatter(range(0, num_regions), num_counts, c=num_counts, marker='o')
            plt.plot(range(0, num_regions), num_counts, "black")
            plt.axis([0, num_regions, 0, max_val])
            plt.ylabel(self.strategy.AXIS_TITLE)
            if ylim:
                plt.ylim(ylim)
        else:
            fig, ax = plt.subplots(len(self.groups[group]),1, sharex=True, squeeze=False)

            max_val = 0
            for i, chrom in enumerate(self.groups[group]):

                num_counts = [self.group_counts[group][track][x] for x in sorted(self.regions) if self.regions[x]["chr"] == chrom]
                num_regions = len(num_counts)

                if max(num_counts) > max_val:
                    max_val = max(num_counts)

                ax[i,0].plot(range(0, num_regions), num_counts, label="g"+str(chrom))
                ax[i,0].text(
                    1.05, 0.5, ("Chr#"+str(chrom)), transform=ax[i,0].transAxes,
                    rotation=270, fontsize=12, va='top',
                    horizontalalignment='center', verticalalignment='center'
                )

                if ylim:
                    ax[i,0].set_ylim(ylim)

            # Y axis label
            fig.text(
                .05, 0.5, self.strategy.AXIS_TITLE, rotation='vertical',
                horizontalalignment='center', verticalalignment='center'
            )

        plt.xlabel("Region# (150bp)")
        plt.suptitle('%s-%s' % (group, track), fontsize=16)

        if annotation:
            plt.annotate(annotation, xy=(.5, 1.03),  xycoords='axes fraction', ha='center', va='center', fontsize=11)

        if save_to:
            plt.savefig(save_to)
            plt.close()
        else:
            plt.show()

    # TODO Copies a lot of plot's functionality
    # TODO Need to alter x ticks to show bin names
    def profile(self, group=None, track="default", ylim=None, save_to=None, annotation=None, bins=None): # pragma: no cover
        """Represent profiled regions in a plot using matplotlib."""

        import matplotlib.pyplot as plt

        #TODO Check this!
        def find_bin(x, bins):
            for i, b in enumerate(sorted(bins)):
                if x < b:
                    if i > 0:
                        return i - 1
                    else:
                        # first bin
                        return 0
            # else last bin
            return len(bins) - 1

        if group is None:
            group = "total"
            fig = plt.subplot(1,1,1)

            if bins:
                num_bins = len(bins)
                bin_contents = np.zeros(len(bins))
                for x in self.group_buckets[group][track]:
                    bin_contents[find_bin(x, bins)] = len(self.group_buckets[group][track][x])
            else:
                print("Unbinned profile not yet supported.")
                import sys
                sys.exit(1)

            max_val = max(bin_contents)
            plt.bar(range(0, num_bins), bin_contents)
            plt.axis([0, num_bins, 0, max_val])
            plt.ylabel("Region Count")


            if ylim:
                plt.ylim(ylim)
        else:
            fig, ax = plt.subplots(len(self.groups[group]),1, sharex=True, squeeze=False)

            for i, chrom in enumerate(self.groups[group]):

                if bins:
                    num_bins = len(bins)
                    bin_contents = np.zeros(len(bins))

                    for x in [self.group_counts[group][track][x] for x in sorted(self.regions) if self.regions[x]["chr"] == chrom]:
                        bin_contents[find_bin(x, bins)] += 1
                else:
                    print("Unbinned profile not yet supported.")
                    import sys
                    sys.exit(1)

                ax[i,0].plot(range(0, num_bins), bin_contents, label="g"+str(chrom))
                ax[i,0].text(
                    1.05, 0.5, ("Chr#"+str(chrom)), transform=ax[i,0].transAxes,
                    rotation=270, fontsize=12, va='top',
                    horizontalalignment='center', verticalalignment='center'
                )

                if ylim:
                    ax[i,0].set_ylim(ylim)

            # Y axis label
            fig.text(
                .05, 0.5, "Region Count", rotation='vertical',
                horizontalalignment='center', verticalalignment='center'
            )

        plt.xlabel("Bin")
        plt.suptitle('%s-%s' % (group, track), fontsize=16)

        if annotation:
            plt.annotate(annotation, xy=(.5, 1.03),  xycoords='axes fraction', ha='center', va='center', fontsize=11)

        if save_to:
            plt.savefig(save_to)
            plt.close()
        else:
            plt.show()

    @property
    def candidates(self):
        if not (len(self.regions) > 0 or len(self.sorted_regions) > 0):
            print("[WARN] No candidates found.\n")

        if self.sorted_regions:
            return [self.regions[i] for i in self.sorted_regions]
        else:
            return self.regions

    #TODO Export to file not stdout!
    #TODO Groupless output
    def export_meta(self, group=None, sep=","):
        """Export census metadata to stdout with the following header: id, track1,
        [track2 ... trackN], chr, pos_start, pos_end. Accepts a seperator but
        defaults to a comma-delimited table."""

        tracks = sorted(self.strategy.TRACKS)

        tracks_headers = []
        groups = sorted(self.groups)
        groups_header = sep.join(groups)
        for g in groups:
            #tracks_headers.append([g + '_' + track for track in tracks])
            tracks_headers.append(sep.join([g + '_' + track for track in tracks]))

        print(sep.join(["id", sep.join(tracks_headers), "chr", "pos_start", "pos_end"]))

        to_iter = self.regions.keys()
        if self.sorted_regions:
            to_iter = self.sorted_regions

        for r in to_iter:
            region = self.regions[r]
            values_string = ""
            for g in groups:
                for t in tracks:
                    values_string += str(self.group_counts[g][t][r])
                    values_string += sep
            values_string = values_string[:-1]
            print(sep.join([
                str(r),
                values_string,
                str(region["chr"]),
                str(region["pos_start"]),
                str(region["pos_end"])]
            ))

    def export_fasta(self, groups=None, filename="out", divide=False):
        """Export all regions held in FASTA format."""

        if not filename.endswith(".fa"):
            filename += ".fa"

        if divide:
            handles = {}
        else:
            current_file = open(filename, "w")

        for region in self:
            if groups is None:
                groups = list(self.__goldilocks.groups.keys())
            for group in groups:
                try:
                    self.__goldilocks.groups[group]
                except KeyError:
                    # Perhaps the user provided a string by mistake...
                    group = groups

                if divide:
                    if group not in handles:
                        handles[group] = open(group + "_" + filename, "w")
                    current_file = handles[group]

                current_file.write(">%s|Chr%s|Pos%d:%d\n" % (group, region["chr"], region["pos_start"], region["pos_end"]))
                current_file.write(self.__goldilocks.groups[group][region["chr"]][region["pos_start"]:region["pos_end"]+1]+"\n")
        if divide:
            for h in handles:
                handles[h].close()
        else:
            current_file.close()


