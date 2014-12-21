import unittest

import numpy as np

from goldilocks.goldilocks import Goldilocks
from goldilocks.strategies import NucleotideCounterStrategy

# TODO Test expected results for max, min, mean, median
class TestGoldilocksRegression_NucleotideCounter(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        sequence_data = {
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
        cls.g = Goldilocks(NucleotideCounterStrategy(["A","C","G","T","N"]), sequence_data, length=3, stride=1)
        cls.EXPECTED_REGIONS = {
                2: {
                    "my_sample": {
                        0: {'A': 1, 'C': 0, 'T': 0, 'G': 0, 'N': 2},
                        1: {'A': 2, 'C': 0, 'T': 0, 'G': 0, 'N': 1},
                        2: {'A': 1, 'C': 0, 'T': 0, 'G': 0, 'N': 2},
                        3: {'A': 2, 'C': 0, 'T': 0, 'G': 0, 'N': 1},
                        4: {'A': 1, 'C': 0, 'T': 0, 'G': 0, 'N': 2},
                        5: {'A': 2, 'C': 0, 'T': 0, 'G': 0, 'N': 1},
                        6: {'A': 1, 'C': 0, 'T': 0, 'G': 0, 'N': 2},
                        7: {'A': 2, 'C': 0, 'T': 0, 'G': 0, 'N': 1},
                    },
                    "my_other_sample": {
                        0: {'A': 1, 'C': 0, 'T': 0, 'G': 1, 'N': 1},
                        1: {'A': 1, 'C': 0, 'T': 0, 'G': 1, 'N': 1},
                        2: {'A': 1, 'C': 0, 'T': 0, 'G': 1, 'N': 1},
                        3: {'A': 1, 'C': 0, 'T': 0, 'G': 1, 'N': 1},
                        4: {'A': 1, 'C': 0, 'T': 0, 'G': 1, 'N': 1},
                        5: {'A': 1, 'C': 0, 'T': 0, 'G': 1, 'N': 1},
                        6: {'A': 1, 'C': 0, 'T': 0, 'G': 1, 'N': 1},
                        7: {'A': 1, 'C': 0, 'T': 0, 'G': 0, 'N': 1},
                    },
                    "total": {
                        0: {'A': 2, 'C': 0, 'T': 0, 'G': 1, 'N': 3, "default": 6},
                        1: {'A': 3, 'C': 0, 'T': 0, 'G': 1, 'N': 2, "default": 6},
                        2: {'A': 2, 'C': 0, 'T': 0, 'G': 1, 'N': 3, "default": 6},
                        3: {'A': 3, 'C': 0, 'T': 0, 'G': 1, 'N': 2, "default": 6},
                        4: {'A': 2, 'C': 0, 'T': 0, 'G': 1, 'N': 3, "default": 6},
                        5: {'A': 3, 'C': 0, 'T': 0, 'G': 1, 'N': 2, "default": 6},
                        6: {'A': 2, 'C': 0, 'T': 0, 'G': 1, 'N': 3, "default": 6},
                        7: {'A': 3, 'C': 0, 'T': 0, 'G': 0, 'N': 2, "default": 5},
                    },
                },
                "X": {
                    "my_sample": {
                        0: {'A': 1, 'C': 0, 'T': 1, 'G': 1, 'N': 0},
                        1: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0},
                        2: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0},
                        3: {'A': 1, 'C': 1, 'T': 1, 'G': 0, 'N': 0},
                        4: {'A': 2, 'C': 1, 'T': 0, 'G': 0, 'N': 0},
                        5: {'A': 1, 'C': 1, 'T': 0, 'G': 1, 'N': 0},
                        6: {'A': 2, 'C': 0, 'T': 0, 'G': 1, 'N': 0},
                        7: {'A': 1, 'C': 0, 'T': 1, 'G': 1, 'N': 0},
                        8: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0},
                        9: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0},
                        10: {'A': 1, 'C': 1, 'T': 1, 'G': 0, 'N': 0},
                        11: {'A': 2, 'C': 1, 'T': 0, 'G': 0, 'N': 0},
                        12: {'A': 1, 'C': 1, 'T': 0, 'G': 0, 'N': 1},
                    },
                    "my_other_sample": {
                        0: {'A': 1, 'C': 0, 'T': 1, 'G': 1, 'N': 0},
                        1: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0},
                        2: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0},
                        3: {'A': 1, 'C': 1, 'T': 1, 'G': 0, 'N': 0},
                        4: {'A': 2, 'C': 1, 'T': 0, 'G': 0, 'N': 0},
                        5: {'A': 1, 'C': 1, 'T': 0, 'G': 1, 'N': 0},
                        6: {'A': 2, 'C': 0, 'T': 0, 'G': 1, 'N': 0},
                        7: {'A': 1, 'C': 0, 'T': 1, 'G': 1, 'N': 0},
                        8: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0},
                        9: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0},
                        10: {'A': 1, 'C': 1, 'T': 1, 'G': 0, 'N': 0},
                        11: {'A': 2, 'C': 1, 'T': 0, 'G': 0, 'N': 0},
                        12: {'A': 1, 'C': 1, 'T': 0, 'G': 0, 'N': 1},
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
                            0: {'A': 1, 'C': 1, 'T': 1, 'G': 0, 'N': 0},
                            1: {'A': 1, 'C': 1, 'T': 1, 'G': 0, 'N': 0},
                            2: {'A': 1, 'C': 1, 'T': 1, 'G': 0, 'N': 0},
                            3: {'A': 1, 'C': 1, 'T': 0, 'G': 0, 'N': 1},
                            4: {'A': 1, 'C': 1, 'T': 0, 'G': 0, 'N': 1},
                            5: {'A': 1, 'C': 1, 'T': 0, 'G': 0, 'N': 1},
                            6: {'A': 1, 'C': 1, 'T': 1, 'G': 0, 'N': 0},
                        },
                        "my_other_sample": {
                            0: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0},
                            1: {'A': 2, 'C': 0, 'T': 1, 'G': 0, 'N': 0},
                            2: {'A': 1, 'C': 0, 'T': 1, 'G': 0, 'N': 1},
                            3: {'A': 1, 'C': 0, 'T': 1, 'G': 0, 'N': 1},
                            4: {'A': 1, 'C': 0, 'T': 1, 'G': 0, 'N': 1},
                            5: {'A': 1, 'C': 0, 'T': 2, 'G': 0, 'N': 0},
                            6: {'A': 2, 'C': 0, 'T': 1, 'G': 0, 'N': 0},
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
                            0: {'A': 1, 'C': 0, 'T': 0, 'G': 0, 'N': 0},
                        },
                        "my_other_sample": {
                            0: {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'N': 1},
                        },
                        "total": {
                            0: {'A': 1, 'C': 0, 'T': 0, 'G': 0, 'N': 1, "default": 2},
                        }
                }
        }

        # 29 regions * 5 bases * (2+1) samples (two samples + total)
        cls.EXPECTED_NUM_REGION = 29
        cls.EXPECTED_REGION_COUNT = cls.EXPECTED_NUM_REGION*5*3

        # Each region gets an additional counter in its total key (total-default)
        cls.EXPECTED_COUNTERS_COUNT = cls.EXPECTED_REGION_COUNT + cls.EXPECTED_NUM_REGION

    def test_group_track_counts_contents(self):
        """Ensure group_counts held by region metadata for each group-track
        combination (including the total-default group-track) match the
        expected number of bases.
        This test is somewhat cheating in that it fetches region metadata from
        the regions dictionary."""

        number_comparisons = 0
        for group in self.g.group_counts:
            for track in self.g.group_counts[group]:
                for region_i, value in enumerate(self.g.group_counts[group][track]):
                    # Get this region_i's chrom and ichr from the region data
                    chrom = self.g.regions[region_i]["chr"]
                    ichr = self.g.regions[region_i]["ichr"]
                    self.assertEqual(self.EXPECTED_REGIONS[chrom][group][ichr][track], value)

                    number_comparisons += 1
                    ichr += 1
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
                        self.assertEqual(self.g.group_counts[group][track][region_id], bucket)
                        number_comparisons += 1
        self.assertEqual(self.EXPECTED_COUNTERS_COUNT, number_comparisons)

class TestGoldilocksRegression_SimpleNucleotideCounter(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.sequence_data = {
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
        cls.g = Goldilocks(NucleotideCounterStrategy(["A","N"]), cls.sequence_data, length=3, stride=3)

        cls.EXPECTED_REGIONS = {
                1: {
                    "my_sample": {
                        0: {'A': 0, 'N': 1},
                        1: {'A': 0, 'N': 1},
                        2: {'A': 0, 'N': 1},
                    },
                    "my_other_sample": {
                        0: {'A': 0, 'N': 1},
                        1: {'A': 0, 'N': 2},
                        2: {'A': 0, 'N': 3},
                    },
                    "total": {
                        0: {'A': 0, 'N': 2, "default": 2},
                        1: {'A': 0, 'N': 3, "default": 3},
                        2: {'A': 0, 'N': 4, "default": 4},
                    },
                },
                2: {
                    "my_sample": {
                        0: {'A': 2, 'N': 0},
                        1: {'A': 2, 'N': 0},
                        2: {'A': 1, 'N': 0},
                    },
                    "my_other_sample": {
                        0: {'A': 1, 'N': 0},
                        1: {'A': 2, 'N': 0},
                        2: {'A': 0, 'N': 0},
                    },
                    "total": {
                        0: {'A': 3, 'N': 0, "default": 3},
                        1: {'A': 4, 'N': 0, "default": 4},
                        2: {'A': 1, 'N': 0, "default": 1},
                    },
                },
                3: {
                    "my_sample": {
                        0: {'A': 0, 'N': 3},
                        1: {'A': 0, 'N': 2},
                        2: {'A': 0, 'N': 0},
                    },
                    "my_other_sample": {
                        0: {'A': 3, 'N': 0},
                        1: {'A': 2, 'N': 0},
                        2: {'A': 0, 'N': 0},
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

        # Each region gets an additional counter in its total key (total-default)
        cls.EXPECTED_COUNTERS_COUNT = cls.EXPECTED_REGION_COUNT + cls.EXPECTED_NUM_REGION

    def test_group_track_counts_contents(self):
        """Ensure group_counts held by region metadata for each group-track
        combination (including the total-default group-track) match the
        expected number of bases.
        This test is somewhat cheating in that it fetches region metadata from
        the regions dictionary."""

        number_comparisons = 0
        for group in self.g.group_counts:
            for track in self.g.group_counts[group]:
                for region_i, value in enumerate(self.g.group_counts[group][track]):
                    # Get this region_i's chrom and ichr from the region data
                    chrom = self.g.regions[region_i]["chr"]
                    ichr = self.g.regions[region_i]["ichr"]
                    self.assertEqual(self.EXPECTED_REGIONS[chrom][group][ichr][track], value)

                    number_comparisons += 1
                    ichr += 1
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
                        self.assertEqual(self.g.group_counts[group][track][region_id], bucket)
                        number_comparisons += 1
        self.assertEqual(self.EXPECTED_COUNTERS_COUNT, number_comparisons)

###############################################################################
# NOTE Tests rely somewhat on chromosomes being processed in numerical order  #
# TODO CandidateLists no longer store values in group_count but do expose g   #
#      which is set to change in future releases... Use 'our' g for tests...  #
###############################################################################
    """
    def test_max_candidates(self):
        EXPECTED_RANK = {
            "A": [4,3,6,7,5,0,1,2,8],
            "N": [2,1,6,0,7,3,4,5,8],
            "default": [6,2,4,7,1,3,0,5,8]
        }
        self.__test_sort_candidates("max", EXPECTED_RANK)

    def test_min_candidates(self):
        EXPECTED_RANK = {
            "A": [0,1,2,8,5,7,3,6,4],
            "N": [3,4,5,8,0,7,1,6,2],
            "default": [8,5,0,1,3,2,4,7,6]
        }
        self.__test_sort_candidates("min", EXPECTED_RANK)
    """

    #TODO Test export_meta for FASTA of each list
    def __test_sort_candidates(self, op, EXPECTED_RANK, targets=None):
        TRACKS = ["A", "N"]

        for group in ["total"]:
            for track in EXPECTED_RANK:
                number_comparisons = 0

                candidates = self.g._filter(op, track=track)
                last_seen = None
                for i, c in enumerate(candidates):
                    # Test value matches expecting region value
                    self.assertEqual(self.g.group_counts[group][track][c["id"]],
                            self.EXPECTED_REGIONS[c["chr"]][group][c["ichr"]][track])

                    # Test region is actually correct
                    total = 0
                    for sample in self.sequence_data:
                        if track == "default":
                            for ttrack in TRACKS:
                                total += self.sequence_data[sample][c["chr"]][c["pos_start"]:c["pos_end"]+1].count(ttrack)
                        else:
                            total += self.sequence_data[sample][c["chr"]][c["pos_start"]:c["pos_end"]+1].count(track)
                    self.assertEqual(self.g.group_counts[group][track][c["id"]], total)

                    # Test expected rank
                    self.assertEqual(EXPECTED_RANK[track][i], c["id"])

                    # Test values are ordered
                    if last_seen is None:
                        last_seen = self.g.group_counts["total"][track][c["id"]]

                    if op == "max":
                        self.assertTrue(self.g.group_counts["total"][track][c["id"]] <= last_seen)
                    elif op == "min":
                        self.assertTrue(self.g.group_counts["total"][track][c["id"]] >= last_seen)
                    elif op == "mean":
                        self.assertTrue(targets[track], candidates._CandidateList__target)
                        if targets is None:
                            self.fail("Invalid test on op:mean using no target.")
                        if track not in targets:
                            self.fail("Invalid test on op:mean using no target.")

                        delta_mean = abs(self.g.group_counts["total"][track][c["id"]] - targets[track])
                        last_delta_mean = abs(last_seen - targets[track])
                        self.assertTrue(delta_mean >= last_delta_mean)
                    else:
                        self.fail("Invalid op.")
                    last_seen = self.g.group_counts["total"][track][c["id"]]

                    number_comparisons += 1

                self.assertEqual(self.EXPECTED_NUM_REGION, number_comparisons)


    def test_mean_candidates(self):
        group = "total"
        EXPECTED_RANK = {}
        EXPECTED_TARGET = {}

        for track in ["A", "N", "default"]:
            total = 0
            count = 0
            for chrom in self.EXPECTED_REGIONS:
                for region in self.EXPECTED_REGIONS[chrom][group]:
                    total += self.EXPECTED_REGIONS[chrom][group][region][track]
                    count += 1
            mean = float(total)/count
            self.assertEqual(self.EXPECTED_NUM_REGION, count)

            count = 0
            delta_mean_buckets = {}
            for chrom in sorted(self.EXPECTED_REGIONS):
                for region in sorted(self.EXPECTED_REGIONS[chrom][group]):
                    delta_mean = abs(self.EXPECTED_REGIONS[chrom][group][region][track] - mean)
                    if delta_mean not in delta_mean_buckets:
                        delta_mean_buckets[delta_mean] = []
                    delta_mean_buckets[delta_mean].append(count)
                    count += 1
            self.assertEqual(self.EXPECTED_NUM_REGION, count)

            argsort_from_mean = []
            for bucket in sorted(delta_mean_buckets):
                argsort_from_mean.extend(delta_mean_buckets[bucket])
            EXPECTED_RANK[track] = argsort_from_mean
            EXPECTED_TARGET[track] = mean

        self.__test_sort_candidates("mean", EXPECTED_RANK, targets=EXPECTED_TARGET)

    def test_median_candidates(self):
        pass

if __name__ == '__main__':
    unittest.main()
