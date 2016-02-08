import unittest

import numpy as np

from goldilocks import Goldilocks
from goldilocks.strategies import NucleotideCounterStrategy
from goldilocks.strategies import GCRatioStrategy

def _setup_sort(suite, op, group, TRACKS):
    EXPECTED_RANK = {}
    EXPECTED_TARGET = {}

    for track in TRACKS:
        total = 0
        count = 0
        scores = []
        for chrom in suite.EXPECTED_REGIONS:
            for region in suite.EXPECTED_REGIONS[chrom][group]:
                scores.append(suite.EXPECTED_REGIONS[chrom][group][region][track])

        if op == "mean":
            target = np.mean(scores)
        elif op == "median":
            target = np.median(scores)
        elif op == "max":
            target = max(scores)
        elif op == "min":
            target = min(scores)
        else:
            suite.fail("Invalid op.")

        suite.assertEqual(suite.EXPECTED_NUM_REGION, len(scores))

        count = 0
        delta_buckets = {}

        # Use sorting in a Python3 friendly fashion...
        chr_num = [chrom for chrom in suite.EXPECTED_REGIONS if type(chrom)==int]
        chr_str = [chrom for chrom in suite.EXPECTED_REGIONS if type(chrom)!=int]
        chroms = sorted(chr_num)
        chroms.extend(sorted(chr_str))

        for chrom in chroms:
            for region in sorted(suite.EXPECTED_REGIONS[chrom][group]):
                delta_target = abs(suite.EXPECTED_REGIONS[chrom][group][region][track] - target)
                if delta_target not in delta_buckets:
                    delta_buckets[delta_target] = []
                delta_buckets[delta_target].append(count)
                count += 1
        suite.assertEqual(suite.EXPECTED_NUM_REGION, count)

        regions_sorted = []
        for bucket in sorted(delta_buckets):
            regions_sorted.extend(delta_buckets[bucket])
        EXPECTED_RANK[track] = regions_sorted
        EXPECTED_TARGET[track] = target

    return EXPECTED_RANK, EXPECTED_TARGET

#TODO Test export_meta for FASTA of each list
def _test_sort_candidates(suite, op, group, track, EXPECTED_RANK, targets=None):
    number_comparisons = 0

    new_g = suite.g.query(op, group=group, track=track)
    candidates = new_g.candidates
    last_seen = None
    print(candidates) # Show some evidence the test is working...
    for i, c in enumerate(candidates):
        # Test value matches expecting region value
        suite.assertEqual(suite.EXPECTED_REGIONS[c["chr"]][group][c["ichr"]][track],
                suite.g.counter_matrix[suite.g._get_group_id(group), suite.g._get_track_id(track), c["id"]])

        # Test region is actually correct
        if suite.g.strategy.RATIO:
            RATIO_OF = suite.g.strategy.RATIO_OF
            if not suite.g.strategy.RATIO_OF:
                RATIO_OF = suite.g.LENGTH
        total = 0
        counter = 0

        if group != "total":
            if track == "default":
                for ttrack in suite.TRACKS:
                    if ttrack == "default" and len(suite.TRACKS) > 1:
                        continue
                    if c["pos_end"] > len(suite.g.groups[group][c["chr"]]):
                        counter += 1
                        continue

                    region = suite.sequence_data[group][c["chr"]][c["pos_start"]-1:c["pos_end"]]
                    if suite.g.strategy.RATIO:
                        total += (suite.g.strategy.census(region, ttrack)) * RATIO_OF
                        counter += 1
                    else:
                        total += suite.g.strategy.census(region, ttrack)
            else:
                if c["pos_end"] > len(suite.g.groups[group][c["chr"]]):
                    counter += 1
                else:
                    region = suite.sequence_data[group][c["chr"]][c["pos_start"]-1:c["pos_end"]]
                    if suite.g.strategy.RATIO:
                        total += (suite.g.strategy.census(region, track)) * RATIO_OF
                        counter += 1
                    else:
                        total += suite.g.strategy.census(region, track)
        else:
            for sample in suite.sequence_data:
                if track == "default":
                    for ttrack in suite.TRACKS:
                        if ttrack == "default" and len(suite.TRACKS) > 1:
                            continue
                        if c["pos_end"] > len(suite.sequence_data[sample][c["chr"]]):
                            counter += 1
                            continue

                        region = suite.sequence_data[sample][c["chr"]][c["pos_start"]-1:c["pos_end"]]
                        if suite.g.strategy.RATIO:
                            total += (suite.g.strategy.census(region, ttrack)) * RATIO_OF
                            counter += 1
                        else:
                            total += suite.g.strategy.census(region, ttrack)
                else:
                    if c["pos_end"] > len(suite.sequence_data[sample][c["chr"]]):
                        counter += 1
                        continue
                    region = suite.sequence_data[sample][c["chr"]][c["pos_start"]-1:c["pos_end"]]
                    if suite.g.strategy.RATIO:
                        total += (suite.g.strategy.census(region, track)) * RATIO_OF
                        counter += 1
                    else:
                        total += suite.g.strategy.census(region, track)

        if suite.g.strategy.RATIO:
            total /= (RATIO_OF * counter)
        suite.assertEqual(suite.g.counter_matrix[suite.g._get_group_id(group), suite.g._get_track_id(track), c["id"]], total)

        # Test expected rank
        suite.assertEqual(EXPECTED_RANK[track][i], c["id"])

        # Test values are ordered
        if last_seen is None:
            last_seen = suite.g.counter_matrix[suite.g._get_group_id(group), suite.g._get_track_id(track), c["id"]]

        if op == "max":
            suite.assertTrue(suite.g.counter_matrix[suite.g._get_group_id(group), suite.g._get_track_id(track), c["id"]] <= last_seen)
        elif op == "min":
            suite.assertTrue(suite.g.counter_matrix[suite.g._get_group_id(group), suite.g._get_track_id(track), c["id"]] >= last_seen)
        elif op == "mean" or op == "median":
            if targets is None:
                suite.fail("Invalid test on op:mean|median using no target.")
            if track not in targets:
                suite.fail("Invalid test on op:mean|median using no target.")

            delta_target = abs(suite.g.counter_matrix[suite.g._get_group_id(group), suite.g._get_track_id(track), c["id"]] - targets[track])
            last_delta_target = abs(last_seen - targets[track])
            suite.assertTrue(delta_target >= last_delta_target)
        else:
            suite.fail("Invalid op.")
        last_seen = suite.g.counter_matrix[suite.g._get_group_id(group), suite.g._get_track_id(track), c["id"]]

        number_comparisons += 1

    suite.assertEqual(suite.EXPECTED_NUM_REGION, number_comparisons)

class TestGoldilocksRegression_NucleotideCounter(unittest.TestCase):

    def reset(self):
        self.sequence_data = {
                "my_sample": {
                    2: "NANANANANA",
                    "X": "GATTACAGATTACAN",
                    "one": "CATCANCAT",
                    "three": "..A",
                },
                "my_other_sample": {
                    2: "GANGANGAN",
                    "X": "GATTACAGATTACAN",
                    "one": "TATANTATA",
                    "three": ".N.",
                }
        }
        self.g = Goldilocks(NucleotideCounterStrategy(["A","C","G","T","N"]), self.sequence_data, length=3, stride=1, ignore_len_mismatch=True)

    def setUp(self):
        self.reset()

    @classmethod
    def setUpClass(cls):
        cls.GROUPS = ["my_sample", "my_other_sample", "total"]
        cls.TRACKS = ["A", "C", "G", "T", "N", "default"]

        cls.EXPECTED_REGIONS = {
                2: {
                    "my_sample": {
                        0: {'A': 1, 'C': 0, 'T': 0, 'G': 0, 'N': 2, "default": 3},
                        1: {'A': 2, 'C': 0, 'T': 0, 'G': 0, 'N': 1, "default": 3},
                        2: {'A': 1, 'C': 0, 'T': 0, 'G': 0, 'N': 2, "default": 3},
                        3: {'A': 2, 'C': 0, 'T': 0, 'G': 0, 'N': 1, "default": 3},
                        4: {'A': 1, 'C': 0, 'T': 0, 'G': 0, 'N': 2, "default": 3},
                        5: {'A': 2, 'C': 0, 'T': 0, 'G': 0, 'N': 1, "default": 3},
                        6: {'A': 1, 'C': 0, 'T': 0, 'G': 0, 'N': 2, "default": 3},
                        7: {'A': 2, 'C': 0, 'T': 0, 'G': 0, 'N': 1, "default": 3},
                    },
                    "my_other_sample": {
                        0: {'A': 1, 'C': 0, 'T': 0, 'G': 1, 'N': 1, "default": 3},
                        1: {'A': 1, 'C': 0, 'T': 0, 'G': 1, 'N': 1, "default": 3},
                        2: {'A': 1, 'C': 0, 'T': 0, 'G': 1, 'N': 1, "default": 3},
                        3: {'A': 1, 'C': 0, 'T': 0, 'G': 1, 'N': 1, "default": 3},
                        4: {'A': 1, 'C': 0, 'T': 0, 'G': 1, 'N': 1, "default": 3},
                        5: {'A': 1, 'C': 0, 'T': 0, 'G': 1, 'N': 1, "default": 3},
                        6: {'A': 1, 'C': 0, 'T': 0, 'G': 1, 'N': 1, "default": 3},
                        7: {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'N': 0, "default": 0},
                    },
                    "total": {
                        0: {'A': 2, 'C': 0, 'T': 0, 'G': 1, 'N': 3, "default": 6},
                        1: {'A': 3, 'C': 0, 'T': 0, 'G': 1, 'N': 2, "default": 6},
                        2: {'A': 2, 'C': 0, 'T': 0, 'G': 1, 'N': 3, "default": 6},
                        3: {'A': 3, 'C': 0, 'T': 0, 'G': 1, 'N': 2, "default": 6},
                        4: {'A': 2, 'C': 0, 'T': 0, 'G': 1, 'N': 3, "default": 6},
                        5: {'A': 3, 'C': 0, 'T': 0, 'G': 1, 'N': 2, "default": 6},
                        6: {'A': 2, 'C': 0, 'T': 0, 'G': 1, 'N': 3, "default": 6},
                        7: {'A': 2, 'C': 0, 'T': 0, 'G': 0, 'N': 1, "default": 3},
                    },
                },
                "X": {
                    "my_sample": {
                        0: {'A': 1, 'C': 0, 'T': 1, 'G': 1, 'N': 0, "default": 3},
                        1: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0, "default": 3},
                        2: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0, "default": 3},
                        3: {'A': 1, 'C': 1, 'T': 1, 'G': 0, 'N': 0, "default": 3},
                        4: {'A': 2, 'C': 1, 'T': 0, 'G': 0, 'N': 0, "default": 3},
                        5: {'A': 1, 'C': 1, 'T': 0, 'G': 1, 'N': 0, "default": 3},
                        6: {'A': 2, 'C': 0, 'T': 0, 'G': 1, 'N': 0, "default": 3},
                        7: {'A': 1, 'C': 0, 'T': 1, 'G': 1, 'N': 0, "default": 3},
                        8: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0, "default": 3},
                        9: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0, "default": 3},
                        10: {'A': 1, 'C': 1, 'T': 1, 'G': 0, 'N': 0, "default": 3},
                        11: {'A': 2, 'C': 1, 'T': 0, 'G': 0, 'N': 0, "default": 3},
                        12: {'A': 1, 'C': 1, 'T': 0, 'G': 0, 'N': 1, "default": 3},
                    },
                    "my_other_sample": {
                        0: {'A': 1, 'C': 0, 'T': 1, 'G': 1, 'N': 0, "default": 3},
                        1: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0, "default": 3},
                        2: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0, "default": 3},
                        3: {'A': 1, 'C': 1, 'T': 1, 'G': 0, 'N': 0, "default": 3},
                        4: {'A': 2, 'C': 1, 'T': 0, 'G': 0, 'N': 0, "default": 3},
                        5: {'A': 1, 'C': 1, 'T': 0, 'G': 1, 'N': 0, "default": 3},
                        6: {'A': 2, 'C': 0, 'T': 0, 'G': 1, 'N': 0, "default": 3},
                        7: {'A': 1, 'C': 0, 'T': 1, 'G': 1, 'N': 0, "default": 3},
                        8: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0, "default": 3},
                        9: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0, "default": 3},
                        10: {'A': 1, 'C': 1, 'T': 1, 'G': 0, 'N': 0, "default": 3},
                        11: {'A': 2, 'C': 1, 'T': 0, 'G': 0, 'N': 0, "default": 3},
                        12: {'A': 1, 'C': 1, 'T': 0, 'G': 0, 'N': 1, "default": 3},
                    },
                    "total": {
                        0: {'A': 2, 'C': 0, 'T': 2, 'G': 2, 'N': 0, "default": 6},
                        1: {'A': 2, 'C': 0, 'T': 4, 'G': 0, 'N': 0, "default": 6},
                        2: {'A': 2, 'C': 0, 'T': 4, 'G': 0, 'N': 0, "default": 6},
                        3: {'A': 2, 'C': 2, 'T': 2, 'G': 0, 'N': 0, "default": 6},
                        4: {'A': 4, 'C': 2, 'T': 0, 'G': 0, 'N': 0, "default": 6},
                        5: {'A': 2, 'C': 2, 'T': 0, 'G': 2, 'N': 0, "default": 6},
                        6: {'A': 4, 'C': 0, 'T': 0, 'G': 2, 'N': 0, "default": 6},
                        7: {'A': 2, 'C': 0, 'T': 2, 'G': 2, 'N': 0, "default": 6},
                        8: {'A': 2, 'C': 0, 'T': 4, 'G': 0, 'N': 0, "default": 6},
                        9: {'A': 2, 'C': 0, 'T': 4, 'G': 0, 'N': 0, "default": 6},
                        10: {'A': 2, 'C': 2, 'T': 2, 'G': 0, 'N': 0, "default": 6},
                        11: {'A': 4, 'C': 2, 'T': 0, 'G': 0, 'N': 0, "default": 6},
                        12: {'A': 2, 'C': 2, 'T': 0, 'G': 0, 'N': 2, "default": 6},
                    },
                },
                "one": {
                        "my_sample": {
                            0: {'A': 1, 'C': 1, 'T': 1, 'G': 0, 'N': 0, "default": 3},
                            1: {'A': 1, 'C': 1, 'T': 1, 'G': 0, 'N': 0, "default": 3},
                            2: {'A': 1, 'C': 1, 'T': 1, 'G': 0, 'N': 0, "default": 3},
                            3: {'A': 1, 'C': 1, 'T': 0, 'G': 0, 'N': 1, "default": 3},
                            4: {'A': 1, 'C': 1, 'T': 0, 'G': 0, 'N': 1, "default": 3},
                            5: {'A': 1, 'C': 1, 'T': 0, 'G': 0, 'N': 1, "default": 3},
                            6: {'A': 1, 'C': 1, 'T': 1, 'G': 0, 'N': 0, "default": 3},
                        },
                        "my_other_sample": {
                            0: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0, "default": 3},
                            1: {'A': 2, 'C': 0, 'T': 1, 'G': 0, 'N': 0, "default": 3},
                            2: {'A': 1, 'C': 0, 'T': 1, 'G': 0, 'N': 1, "default": 3},
                            3: {'A': 1, 'C': 0, 'T': 1, 'G': 0, 'N': 1, "default": 3},
                            4: {'A': 1, 'C': 0, 'T': 1, 'G': 0, 'N': 1, "default": 3},
                            5: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0, "default": 3},
                            6: {'A': 2, 'C': 0, 'T': 1, 'G': 0, 'N': 0, "default": 3},
                        },
                        "total": {
                            0: {'A': 2, 'C': 1, 'T': 3, 'G': 0, 'N': 0, "default": 6},
                            1: {'A': 3, 'C': 1, 'T': 2, 'G': 0, 'N': 0, "default": 6},
                            2: {'A': 2, 'C': 1, 'T': 2, 'G': 0, 'N': 1, "default": 6},
                            3: {'A': 2, 'C': 1, 'T': 1, 'G': 0, 'N': 2, "default": 6},
                            4: {'A': 2, 'C': 1, 'T': 1, 'G': 0, 'N': 2, "default": 6},
                            5: {'A': 2, 'C': 1, 'T': 2, 'G': 0, 'N': 1, "default": 6},
                            6: {'A': 3, 'C': 1, 'T': 2, 'G': 0, 'N': 0, "default": 6},
                        },
                },
                "three": {
                        "my_sample": {
                            0: {'A': 1, 'C': 0, 'T': 0, 'G': 0, 'N': 0, "default": 1},
                        },
                        "my_other_sample": {
                            0: {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'N': 1, "default": 1},
                        },
                        "total": {
                            0: {'A': 1, 'C': 0, 'T': 0, 'G': 0, 'N': 1, "default": 2},
                        }
                }
        }

        # 29 regions * 5 bases * (2+1) samples (two samples + total)
        cls.EXPECTED_NUM_REGION = 29
        cls.EXPECTED_REGION_COUNT = cls.EXPECTED_NUM_REGION*5*3

        # Each region gets an additional counter
        cls.EXPECTED_COUNTERS_COUNT = cls.EXPECTED_REGION_COUNT + cls.EXPECTED_NUM_REGION*3

    def test_group_track_counts_contents(self):
        """Ensure group_counts held by region metadata for each group-track
        combination (including the total-default group-track) match the
        expected number of bases.
        This test is somewhat cheating in that it fetches region metadata from
        the regions dictionary."""

        number_comparisons = 0
        ggroups = self.g.groups.keys() + ["total"]
        ttracks = self.g.strategy.TRACKS + ["default"]
        for group in ggroups:
            for track in ttracks:
                for region_i, value in enumerate(self.g.counter_matrix[self.g._get_group_id(group), self.g._get_track_id(track),]):
                    # Get this region_i's chrom and ichr from the region data
                    chrom = self.g.regions[region_i]["chr"]
                    ichr = self.g.regions[region_i]["ichr"]
                    self.assertEqual(self.EXPECTED_REGIONS[chrom][group][ichr][track], value)

                    number_comparisons += 1
        self.assertEqual(self.EXPECTED_COUNTERS_COUNT, number_comparisons)

    def test_group_track_bucket_contents(self):
        """Check that regions appear in the correct group_buckets for each
        group-track combination (including the total-default group-track).
        The test is somewhat of a cheat in that it assumes the contents of the
        group_counts counters are correct as tested in
        test_group_track_counts_contents."""

        number_comparisons = 0
        for group in self.g.group_buckets:
            for track in self.g.group_buckets[group]:
                for bucket in self.g.group_buckets[group][track]:
                    for region_id in self.g.group_buckets[group][track][bucket]:
                        self.assertEqual(self.g.counter_matrix[self.g._get_group_id(group), self.g._get_track_id(track), region_id], bucket)
                        number_comparisons += 1
        self.assertEqual(self.EXPECTED_COUNTERS_COUNT, number_comparisons)

    def test_max_candidates(self):
        for group in self.GROUPS:
            for track in self.TRACKS:
                self.reset()
                EXPECTED_RANK, EXPECTED_TARGET = _setup_sort(self, "max", group, self.TRACKS)
                _test_sort_candidates(self, "max", group, track, EXPECTED_RANK)

    def test_min_candidates(self):
        for group in self.GROUPS:
            for track in self.TRACKS:
                self.reset()
                EXPECTED_RANK, EXPECTED_TARGET = _setup_sort(self, "min", group, self.TRACKS)
                _test_sort_candidates(self, "min", group, track, EXPECTED_RANK)

    def test_mean_candidates(self):
        for group in self.GROUPS:
            for track in self.TRACKS:
                self.reset()
                EXPECTED_RANK, EXPECTED_TARGET = _setup_sort(self, "mean", group, self.TRACKS)
                _test_sort_candidates(self, "mean", group, track, EXPECTED_RANK, targets=EXPECTED_TARGET)

    def test_median_candidates(self):
        for group in self.GROUPS:
            for track in self.TRACKS:
                self.reset()
                EXPECTED_RANK, EXPECTED_TARGET = _setup_sort(self, "median", group, self.TRACKS)
                _test_sort_candidates(self, "median", group, track, EXPECTED_RANK, targets=EXPECTED_TARGET)

#TODO Test export_meta
#TODO Test percentile_distance around ops
#TODO Test actual_distance around ops
class TestGoldilocksRegression_SimpleNucleotideCounter(unittest.TestCase):
    def setUp(self):
        self.reset()

    def reset(self):
        self.sequence_data = {
                "my_sample": {
                    1: "..N..N..N",
                    2: "A.A.AA..A",
                    3: "NNN.NN...",
                },
                "my_other_sample": {
                    1: "N..NN.NNN",
                    2: "A..AA....",
                    3: "AAA.AA...",
                }
        }
        self.g = Goldilocks(NucleotideCounterStrategy(["A","N"]), self.sequence_data, length=3, stride=3)

    @classmethod
    def setUpClass(cls):
        cls.GROUPS = ["my_sample", "my_other_sample", "total"]
        cls.TRACKS = ["A", "N", "default"]

        cls.EXPECTED_REGIONS = {
                1: {
                    "my_sample": {
                        0: {'A': 0, 'N': 1, "default": 1},
                        1: {'A': 0, 'N': 1, "default": 1},
                        2: {'A': 0, 'N': 1, "default": 1},
                    },
                    "my_other_sample": {
                        0: {'A': 0, 'N': 1, "default": 1},
                        1: {'A': 0, 'N': 2, "default": 2},
                        2: {'A': 0, 'N': 3, "default": 3},
                    },
                    "total": {
                        0: {'A': 0, 'N': 2, "default": 2},
                        1: {'A': 0, 'N': 3, "default": 3},
                        2: {'A': 0, 'N': 4, "default": 4},
                    },
                },
                2: {
                    "my_sample": {
                        0: {'A': 2, 'N': 0, "default": 2},
                        1: {'A': 2, 'N': 0, "default": 2},
                        2: {'A': 1, 'N': 0, "default": 1},
                    },
                    "my_other_sample": {
                        0: {'A': 1, 'N': 0, "default": 1},
                        1: {'A': 2, 'N': 0, "default": 2},
                        2: {'A': 0, 'N': 0, "default": 0},
                    },
                    "total": {
                        0: {'A': 3, 'N': 0, "default": 3},
                        1: {'A': 4, 'N': 0, "default": 4},
                        2: {'A': 1, 'N': 0, "default": 1},
                    },
                },
                3: {
                    "my_sample": {
                        0: {'A': 0, 'N': 3, "default": 3},
                        1: {'A': 0, 'N': 2, "default": 2},
                        2: {'A': 0, 'N': 0, "default": 0},
                    },
                    "my_other_sample": {
                        0: {'A': 3, 'N': 0, "default": 3},
                        1: {'A': 2, 'N': 0, "default": 2},
                        2: {'A': 0, 'N': 0, "default": 0},
                    },
                    "total": {
                        0: {'A': 3, 'N': 3, "default": 6},
                        1: {'A': 2, 'N': 2, "default": 4},
                        2: {'A': 0, 'N': 0, "default": 0},
                    },
                }
        }

        # 9 regions * 2 bases * (2+1) samples (two samples + total)
        cls.EXPECTED_NUM_REGION = 9
        cls.EXPECTED_REGION_COUNT = cls.EXPECTED_NUM_REGION*2*3

        # Each region gets an additional default counter
        cls.EXPECTED_COUNTERS_COUNT = cls.EXPECTED_REGION_COUNT + cls.EXPECTED_NUM_REGION*3

    def test_group_track_counts_contents(self):
        """Ensure group_counts held by region metadata for each group-track
        combination (including the total-default group-track) match the
        expected number of bases.
        This test is somewhat cheating in that it fetches region metadata from
        the regions dictionary."""

        number_comparisons = 0
        ggroups = self.g.groups.keys() + ["total"]
        ttracks = self.g.strategy.TRACKS + ["default"]
        for group in ggroups:
            for track in ttracks:
                for region_i, value in enumerate(self.g.counter_matrix[self.g._get_group_id(group), self.g._get_track_id(track),]):
                    # Get this region_i's chrom and ichr from the region data
                    chrom = self.g.regions[region_i]["chr"]
                    ichr = self.g.regions[region_i]["ichr"]
                    self.assertEqual(self.EXPECTED_REGIONS[chrom][group][ichr][track], value)

                    number_comparisons += 1
        self.assertEqual(self.EXPECTED_COUNTERS_COUNT, number_comparisons)

    def test_group_track_bucket_contents(self):
        """Check that regions appear in the correct group_buckets for each
        group-track combination (including the total-default group-track).
        The test is somewhat of a cheat in that it assumes the contents of the
        group_counts counters are correct as tested in
        test_group_track_counts_contents."""

        number_comparisons = 0
        for group in self.g.group_buckets:
            for track in self.g.group_buckets[group]:
                for bucket in self.g.group_buckets[group][track]:
                    for region_id in self.g.group_buckets[group][track][bucket]:
                        self.assertEqual(self.g.counter_matrix[self.g._get_group_id(group), self.g._get_track_id(track), region_id], bucket)
                        number_comparisons += 1
        self.assertEqual(self.EXPECTED_COUNTERS_COUNT, number_comparisons)

    def test_max_candidates(self):
        for group in self.GROUPS:
            for track in self.TRACKS:
                self.reset()
                EXPECTED_RANK, EXPECTED_TARGET = _setup_sort(self, "max", group, self.TRACKS)
                _test_sort_candidates(self, "max", group, track, EXPECTED_RANK)

    def test_min_candidates(self):
        for group in self.GROUPS:
            for track in self.TRACKS:
                self.reset()
                EXPECTED_RANK, EXPECTED_TARGET = _setup_sort(self, "min", group, self.TRACKS)
                _test_sort_candidates(self, "min", group, track, EXPECTED_RANK)

    def test_mean_candidates(self):
        for group in self.GROUPS:
            for track in self.TRACKS:
                self.reset()
                EXPECTED_RANK, EXPECTED_TARGET = _setup_sort(self, "mean", group, self.TRACKS)
                _test_sort_candidates(self, "mean", group, track, EXPECTED_RANK, targets=EXPECTED_TARGET)

    def test_median_candidates(self):
        for group in self.GROUPS:
            for track in self.TRACKS:
                self.reset()
                EXPECTED_RANK, EXPECTED_TARGET = _setup_sort(self, "median", group, self.TRACKS)
                _test_sort_candidates(self, "median", group, track, EXPECTED_RANK, targets=EXPECTED_TARGET)


class TestGoldilocksRegression_SimpleGCRatioCounter(unittest.TestCase):

    def setUp(self):
        self.reset()

    def reset(self):
        self.sequence_data = {
                "my_sample": {
                    1: "GCGCGCGC..GCGCGC....GCGC......GC",
                },
                "my_other_sample": {
                    1: "GC......GCGC....GCGCGC..GCGCGCGCAAAAAAAA",
                }
        }
        self.g = Goldilocks(GCRatioStrategy(), self.sequence_data, length=8, stride=8, ignore_len_mismatch=True)

    @classmethod
    def setUpClass(cls):
        cls.GROUPS = ["my_sample", "my_other_sample", "total"]
        cls.TRACKS = ["default"]

        cls.EXPECTED_REGIONS = {
                1: {
                    "my_sample": {
                        0: {"default": 1.0},
                        1: {"default": 0.75},
                        2: {"default": 0.5},
                        3: {"default": 0.25},
                        4: {"default": 0},
                    },
                    "my_other_sample": {
                        0: {"default": 0.25},
                        1: {"default": 0.5},
                        2: {"default": 0.75},
                        3: {"default": 1.0},
                        4: {"default": 0},
                    },
                    "total": {
                        0: {"default": 0.625},
                        1: {"default": 0.625},
                        2: {"default": 0.625},
                        3: {"default": 0.625},
                        4: {"default": 0},
                    },
                },
        }

        cls.EXPECTED_NUM_REGION = 5
        cls.EXPECTED_REGION_COUNT = cls.EXPECTED_NUM_REGION*3

        # Each region gets an additional default counter
        cls.EXPECTED_COUNTERS_COUNT = cls.EXPECTED_REGION_COUNT

    def test_group_track_counts_contents(self):
        """Ensure group_counts held by region metadata for each group-track
        combination (including the total-default group-track) match the
        expected number of bases.
        This test is somewhat cheating in that it fetches region metadata from
        the regions dictionary."""

        number_comparisons = 0
        ggroups = self.g.groups.keys() + ["total"]
        ttracks = self.g.strategy.TRACKS
        for group in ggroups:
            for track in ttracks:
                for region_i, value in enumerate(self.g.counter_matrix[self.g._get_group_id(group), self.g._get_track_id(track),]):
                    # Get this region_i's chrom and ichr from the region data
                    chrom = self.g.regions[region_i]["chr"]
                    ichr = self.g.regions[region_i]["ichr"]
                    self.assertEqual(self.EXPECTED_REGIONS[chrom][group][ichr][track], value)

                    number_comparisons += 1
        self.assertEqual(self.EXPECTED_COUNTERS_COUNT, number_comparisons)

    def test_group_track_bucket_contents(self):
        """Check that regions appear in the correct group_buckets for each
        group-track combination (including the total-default group-track).
        The test is somewhat of a cheat in that it assumes the contents of the
        group_counts counters are correct as tested in
        test_group_track_counts_contents."""

        number_comparisons = 0
        for group in self.g.group_buckets:
            for track in self.g.group_buckets[group]:
                for bucket in self.g.group_buckets[group][track]:
                    for region_id in self.g.group_buckets[group][track][bucket]:
                        self.assertEqual(self.g.counter_matrix[self.g._get_group_id(group), self.g._get_track_id(track), region_id], bucket)
                        number_comparisons += 1
        self.assertEqual(self.EXPECTED_COUNTERS_COUNT, number_comparisons)

    def test_max_candidates(self):
        for group in self.GROUPS:
            for track in self.TRACKS:
                self.reset()
                EXPECTED_RANK, EXPECTED_TARGET = _setup_sort(self, "max", group, self.TRACKS)
                _test_sort_candidates(self, "max", group, track, EXPECTED_RANK)

    def test_min_candidates(self):
        for group in self.GROUPS:
            for track in self.TRACKS:
                self.reset()
                EXPECTED_RANK, EXPECTED_TARGET = _setup_sort(self, "min", group, self.TRACKS)
                _test_sort_candidates(self, "min", group, track, EXPECTED_RANK)

    def test_mean_candidates(self):
        for group in self.GROUPS:
            for track in self.TRACKS:
                self.reset()
                EXPECTED_RANK, EXPECTED_TARGET = _setup_sort(self, "mean", group, self.TRACKS)
                _test_sort_candidates(self, "mean", group, track, EXPECTED_RANK, targets=EXPECTED_TARGET)

    def test_median_candidates(self):
        for group in self.GROUPS:
            for track in self.TRACKS:
                self.reset()
                EXPECTED_RANK, EXPECTED_TARGET = _setup_sort(self, "median", group, self.TRACKS)
                _test_sort_candidates(self, "median", group, track, EXPECTED_RANK, targets=EXPECTED_TARGET)

if __name__ == '__main__':
    unittest.main()
