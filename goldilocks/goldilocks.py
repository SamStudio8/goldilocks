__author__ = "Sam Nicholls <sn8@sanger.ac.uk>"
__copyright__ = "Copyright (c) Sam Nicholls"
__version__ = "0.0.5"
__maintainer__ = "Sam Nicholls <sam@samnicholls.net>"

import numpy as np
from math import floor, ceil

class CandidateList(list):
    """A list defining its own tab-delimited table string representation when
    printed by a user. Provides other utility methods for generating other useful
    outputs including exporting sequences to FASTA."""
    def __init__(self, g, group, *args):
        self.__goldilocks = g
        if group is None:
            self.__group = "total"
        else:
            self.__group = group
        list.__init__(self, *args)

    #TODO I don't like how we're keeping a reference to the Goldilocks object
    #     here, in future I'll add a method to populate some form of DataFrame
    #     that contains all the data required without needing Goldilocks behind.
    def __repr__(self):
        str_rep = "ID\tVAL\tCHR\tPOSITIONS (INC.)\n"
        for region in self:
            str_rep += ("%d\t%s\t%s\t%10d - %10d\n" % (region["id"],
                                            self.__goldilocks.group_counts[self.__group]["default"][region["id"]],
                                            region["chr"],
                                            region["pos_start"],
                                            region["pos_end"],
            ))
        return str_rep

    #TODO Export to file not stdout!
    #TODO Support groupless requests - export each group to seperate file or all together
    def export_fasta(self, group):
        """Export all regions held in the CandidateList in FASTA format."""
        for region in self:
            print(">%s|Chr%s|Pos%d:%d" % (group, region["chr"], region["pos_start"], region["pos_end"]))
            print(self.__goldilocks.groups[group][region["chr"]][region["pos_start"]:region["pos_end"]+1])


# TODO Generate database of regions with stats... SQL/SQLite
#      - Probably more of a wrapper script than core-functionality: goldib
class Goldilocks(object):
    """Facade class responsible for conducting a census of provided genomic regions
    using the given strategy and provides an interface via _filter to query results
    for given criteria and return a CandidateList."""

    def __load_chromosome(self, arr, data, track):
        return self.strategy.prepare(arr, data, track)

    def __init__(self, strategy, data, is_seq=True, length=None, stride=1, med_window=12.5):

        self.strategy = strategy# Search strategy

        self.chr_max_len = {}   # Map chromosomes to the largest variant position
                                # seen across all files

        self.groups = {}        # For each group stores a dict with chromosome
                                # numbers as keys with lists of variant positions
                                # as values

        self.group_buckets = {} # For each group holds keys of region sizes with
                                # values a list of region_i of that size

        self.group_counts = {}  # Holds a list of each region size seen for each
                                # group for calculating quantiles etc.

        self.regions = {}       # Stores information on each region checked

        self.MED_WINDOW = med_window # Middle 25%, can also be overriden later

        # Read data
        self.groups = data
        self.max_chr_max_len = None

        for group in self.groups:
            #TODO Catch duplicates etc..
            self.group_buckets[group] = {}
            self.group_counts[group] = {}

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

                if self.max_chr_max_len is None:
                    self.max_chr_max_len = len_current_seq
                elif len_current_seq > self.max_chr_max_len:
                    self.max_chr_max_len = len_current_seq

        self.group_buckets["total"] = {}
        self.group_counts["total"] = {}

        # Ensure stride and length are valid, is a length has not been provided
        # use the size of the largest chromosome divided by 10.
        if stride < 1:
            raise Exception("[FAIL] Stride must be at least 1 base wide.")
        self.STRIDE = stride

        # TODO Somewhat arbitrary...
        if length is None:
            length = self.max_chr_max_len / 10

        if length < 1:
            raise Exception("[FAIL] Length must be at least 1 base wide.")
        self.LENGTH = length

        # Conduct a census of the regions on each chromosome using the user
        # defined length and stride. Counting the number of variants present for
        # each group.
        regions = {}
        region_i = 0
        for chrno, size in sorted(self.chr_max_len.items()):
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
                    "ichr": i,
                    "chr": chrno,
                    "pos_start": region_s,
                    "pos_end": region_e
                }

                for group in self.groups:
                    for track in self.strategy.TRACKS:
                        value = self.strategy.evaluate(chros[group][track][region_s:region_e+1], track=track)

                        # TODO Should we be ignoring these regions if they are empty?
                        # TODO Config option to include 0 in flter metrics
                        if track not in self.group_buckets[group]:
                            self.group_buckets[group][track] = {}

                        if value not in self.group_buckets[group][track]:
                            # Add this particular number of variants as a bucket
                            self.group_buckets[group][track][value] = []

                        # Add the region id to the bucket
                        self.group_buckets[group][track][value].append(region_i)

#                   if value > 0:
                        # Append the number of variants counted in this region
                        # for this group to a list
                        if track not in self.group_counts["total"]:
                            self.group_counts["total"][track] = []
                        if track not in self.group_counts[group]:
                            self.group_counts[group][track] = []

                        self.group_counts[group][track].append(value)

                        try:
                            self.group_counts["total"][track][region_i] += value
                        except IndexError:
                            self.group_counts["total"][track].append(value)

                        # Initialise better...
                        if (len(self.strategy.TRACKS) > 1
                                or (len(self.strategy.TRACKS) == 1
                                    and "default" not in self.strategy.TRACKS)):
                            #TODO Halt when creating a strategy track called default
                            if "default" not in self.group_counts["total"]:
                                self.group_counts["total"]["default"] = []

                            try:
                                self.group_counts["total"]["default"][region_i] += value
                            except IndexError:
                                self.group_counts["total"]["default"].append(value)

                region_i += 1

        # Populate total group_buckets now census is complete
        self.group_buckets["total"] = {}
        super_totals = []
        for track in self.strategy.TRACKS:
            totals = []
            for region_id, value in enumerate(self.group_counts["total"][track]):
                try:
                    totals[region_id] += value
                except IndexError:
                    totals.append(value)

                try:
                    super_totals[region_id] += value
                except IndexError:
                    super_totals.append(value)

            total_buckets = {}
            for region_id, total in enumerate(totals):
                if total not in total_buckets:
                    total_buckets[total] = []
                total_buckets[total].append(region_id)

            self.group_buckets["total"][track] = total_buckets

        # Populate total-default group-track which sums the totals across
        # all tracks...
        total_buckets = {}
        for region_id, total in enumerate(super_totals):
            if total not in total_buckets:
                total_buckets[total] = []
            total_buckets[total].append(region_id)
        self.group_buckets["total"]["default"] = total_buckets

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
            if actual:
                q_low  = np.percentile(np.asarray(self.group_counts[group][track]), 50) - lower_window
                q_high = np.percentile(np.asarray(self.group_counts[group][track]), 50) + upper_window
            else:
                q_low  = np.percentile(np.asarray(self.group_counts[group][track]), 50 - lower_window)
                q_high = np.percentile(np.asarray(self.group_counts[group][track]), 50 + upper_window)
            target = np.percentile(np.asarray(self.group_counts[group][track]), 50)
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
            raise NotImplementedError("[FAIL] Function '%s' not supported" % func.lower())
        return float(q_low), float(q_high), float(target)


    def __check_exclusions(self, exclusions, region_dict, use_and=False, use_chrom=False):
        if exclusions is None or len(exclusions) == 0:
            return False

        def exclude_chro(region_dict, chr_list_or_bool):
            try:
                if region_dict["chr"] in chr_list_or_bool:
                    return True
            except TypeError:
                # It's probably a bool not a list (it could be an int if users
                # have failed to read the documentation...)
                if chr_list_or_bool is True:
                    return True
            return False

        def exclude_start(region_dict, operand, position):
            if operand < 0:
                if region_dict["pos_start"] <= position:
                    return True
            elif operand > 0:
                if region_dict["pos_start"] >= position:
                    return True

            return False

        def exclude_end(region_dict, operand, position):
            if operand < 0:
                if region_dict["pos_end"] <= position:
                    return True
            elif operand > 0:
                if region_dict["pos_end"] >= position:
                    return True

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
                ret = exclude_chro(region_dict, to_apply["chr"])
            elif name == "start_lte":
                ret = exclude_start(region_dict, -1, to_apply["start_lte"])
            elif name == "start_gte":
                ret = exclude_start(region_dict, 1, to_apply["start_gte"])
            elif name == "end_lte":
                ret = exclude_end(region_dict, -1, to_apply["end_lte"])
            elif name == "end_gte":
                ret = exclude_end(region_dict, 1, to_apply["end_gte"])
            else:
                #TODO Better handling of invalid exclusion property
                print("[WARN] Attempted to exclude on invalid property '%s'" % name)

            if use_and:
                # Require all exclusions to be true...
                if ret is False:
                    return False
            else:
                # If we're not waiting on all conditions, we can exclude on the first
                if ret is True:
                    return True

        if use_and:
            # If we didn't bail on a previous false, all conditions must be satisfied
            return True
        return False

    # TODO Pretty hacky at the moment... Just trying some things out!
    def _filter(self, func="median", track="default", actual_distance=None, percentile_distance=None,
            direction=0, group=None, limit=0, exclusions=None, use_and=False, use_chrom=False):
        distance = None
        actual = False
        if actual_distance is not None and percentile_distance is not None:
            raise Exception("[FAIL] Cannot filter by both actual_distance and percentile_difference. Select one.")

        if actual_distance is not None:
            distance = actual_distance
            actual = True

        if percentile_distance is not None:
            distance = percentile_distance

        if group is None:
            group = "total"

        candidates = []
        if distance is not None:
            upper_window = float(distance)
            lower_window = float(distance)
            if func not in ["max", "min"]:
                upper_window = float(distance) / 2
                lower_window = float(distance) / 2

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

        filtered = CandidateList(self, group)
        for region in sorted(self.regions,
                        key=lambda x: (abs(self.group_counts[group][track][x] - target),
                            self.group_counts[group][track][x])):
            num_total += 1
            if region in candidates:
                num_selected += 1
                if not self.__check_exclusions(exclusions, self.regions[region], use_and, use_chrom):
                    self.regions[region]["id"] = region
                    filtered.append(self.regions[region])
                else:
                    num_excluded += 1

        # Return the top N elements if desired
        # TODO Report total, selected, selected-excluded and selected-filtered
        if limit:
            filtered = CandidateList(self, group, filtered[0:limit])

        print("[NOTE] %d processed, %d match search criteria, %d excluded, %d limit" %
                (num_total, num_selected, num_excluded, limit))

        return filtered

    def plot(self, group=None, track="default", ylim=None, save_to=None, annotation=None):
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

    #TODO Export to file not stdout!
    #TODO Groupless output
    def export_meta(self, group, sep=","):
        """Export census metadata to stdout with the following header: id, track1,
        [track2 ... trackN], chr, pos_start, pos_end. Accepts a seperator but
        defaults to a comma-delimited table."""

        tracks = sorted(self.strategy.TRACKS)
        tracks_header = sep.join(tracks)

        print(sep.join(["id", tracks_header, "chr", "pos_start", "pos_end"]))
        for r in self.regions:
            region = self.regions[r]
            values_string = ""
            for t in tracks:
                values_string += str(self.group_counts[group][t][r])
                values_string += sep
            values_string = values_string[:-1]
            print(sep.join([
                str(r),
                values_string,
                str(region["chr"]),
                str(region["pos_start"]),
                str(region["pos_end"])]
            ))
