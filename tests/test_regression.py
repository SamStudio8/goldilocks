import unittest

import numpy as np

from goldilocks import Goldilocks
from goldilocks.strategies import NucleotideCounterStrategy

################################################################################
# NOTE Tests following are hard coded regression tests whose answers were      #
#      generated and have been checked 'by eye'. Such tests are designed to    #
#      ensure that regressions have not been introduced in to Goldilocks,      #
#      rather than to verify accuracy of its results.                          #
################################################################################

DATA = {
    "GroupOne": {
        1: "AAACCCNNNTTTAAACCCGGGNNNAAACCCGGGTTTNNNCCCGGGTTT",
        2: "GCGTNANNGANGGCTANTCTANAGCNNTTTCTNTNNGCANCANTTGNN",
    },
    "GroupTwo": {
        1: "AAACCCNNNTTTAAACCCGGGNNNAAACCCGGGTTTNNNCCCGGGTTT",
        2: "GCGTNANNGANGGCTANTCTANAGCNNTTTCTNTNNGCANCANTTGNN",
    }
}
DATA_FAI = {
    "GroupOne": {
        "file": "tests/dat/my_sequence.fa.fai"
    },
    "GroupTwo": {
        "file": "tests/dat/my_other_sequence.fa.fai"
    }
}
STRIDE = 1
LENGTH = 3

# TODO Test group_counts
# TODO Hash Graphs and FASTA
# TODO Test sorting and filtering beyond min
# TODO Test regions for GroupTwo/total
class TestGoldilocksRegression_NCounter(unittest.TestCase):

#  GroupOne
#  CHR 1
#
#       1 |======***============***============***=========| 48
#          AAACCCNNNTTTAAACCCGGGNNNAAACCCGGGTTTNNNCCCGGGTTT
#          |0| |1| |1| |0| |0| |2| |0| |0| |0| |3| |0| |0|
#           |0| |2| |0| |0| |0| |3| |0| |0| |0| |2| |0| |0|
#            |0| |3| |0| |0| |0| |2| |0| |0| |1| |1| |0|
#             |0| |2| |0| |0| |1| |1| |0| |0| |2| |0| |0|
#

#  GroupOne
#  CHR 2
#
#       1 |====*=**==*=====*====*===**=====*=**===*==*===**| 48
#          GCGTNANNGANGGCTANTCTANAGCNNTTTCTNTNNGCANCANTTGNN
#          |0| |2| |1| |0| |1| |1| |2| |0| |2| |0| |1| |1|
#           |0| |2| |1| |0| |0| |1| |2| |0| |2| |1| |1| |2|
#            |1| |2| |1| |1| |0| |0| |1| |1| |2| |1| |1|
#             |1| |1| |0| |1| |1| |1| |0| |1| |1| |1| |0|
#
    @classmethod
    def setUpClass(cls):
        cls.g = Goldilocks(NucleotideCounterStrategy(["N"]), DATA, length=LENGTH, stride=STRIDE)
        cls.EXPECTED_REGIONS = {
            1: {
                0: 0,
                1: 0,
                2: 0,
                3: 0,
                4: 1,
                5: 2,
                6: 3,
                7: 2,
                8: 1,
                9: 0,
                10: 0,
                11: 0,
                12: 0,
                13: 0,
                14: 0,
                15: 0,
                16: 0,
                17: 0,
                18: 0,
                19: 1,
                20: 2,
                21: 3,
                22: 2,
                23: 1,
                24: 0,
                25: 0,
                26: 0,
                27: 0,
                28: 0,
                29: 0,
                30: 0,
                31: 0,
                32: 0,
                33: 0,
                34: 1,
                35: 2,
                36: 3,
                37: 2,
                38: 1,
                39: 0,
                40: 0,
                41: 0,
                42: 0,
                43: 0,
                44: 0,
                45: 0,
            },
            2: {
                0: 0,
                1: 0,
                2: 1,
                3: 1,
                4: 2,
                5: 2,
                6: 2,
                7: 1,
                8: 1,
                9: 1,
                10: 1,
                11: 0,
                12: 0,
                13: 0,
                14: 1,
                15: 1,
                16: 1,
                17: 0,
                18: 0,
                19: 1,
                20: 1,
                21: 1,
                22: 0,
                23: 1,
                24: 2,
                25: 2,
                26: 1,
                27: 0,
                28: 0,
                29: 0,
                30: 1,
                31: 1,
                32: 2,
                33: 2,
                34: 2,
                35: 1,
                36: 0,
                37: 1,
                38: 1,
                39: 1,
                40: 1,
                41: 1,
                42: 1,
                43: 0,
                44: 1,
                45: 2,
            }
        }

        cls.EXPECTED_REGION_TOTALS = []
        for j in cls.EXPECTED_REGIONS:
            for i in cls.EXPECTED_REGIONS[j]:
                cls.EXPECTED_REGION_TOTALS.append(cls.EXPECTED_REGIONS[j][i]*2)

        cls.EXPECTED_BUCKETS = {
            0: [0, 1, 2, 3, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 24, 25, 26,
                    27, 28, 29, 30, 31, 32, 33, 39, 40, 41, 42, 43, 44, 45,
                        46, 47, 57, 58, 59, 63, 64, 68, 73, 74, 75, 82, 89],
            1: [4, 8, 19, 23, 34, 38,
                        48, 49, 53, 54, 55, 56, 60, 61, 62, 65, 66, 67, 69,
                        72, 76, 77, 81, 83, 84, 85, 86, 87, 88, 90],
            2: [5, 7, 20, 22, 35, 37,
                        50, 51, 52, 70, 71, 78, 79, 80, 91],
            3: [6, 21, 36],
        }
        cls.EXPECTED_REGION_COUNT = len(cls.EXPECTED_REGIONS[1]) + len(cls.EXPECTED_REGIONS[2])

    def test_number_regions(self):
        self.assertEqual(self.EXPECTED_REGION_COUNT, len(self.g.regions.items()))

    # TODO Technically these check the metadata is correct and not the *actual* regions
    def test_region_lengths(self):
        number_comparisons = 0
        # Ensure ALL meet LENGTH
        for region_i, region_data in self.g.regions.items():
            number_comparisons += 1
            self.assertEqual(LENGTH, len(range(region_data["pos_start"], region_data["pos_end"])) + 1)
        self.assertEqual(self.EXPECTED_REGION_COUNT, number_comparisons)

    # TODO Technically these check the metadata is correct and not the *actual* regions
    def test_region_stride(self):
        number_comparisons = 0
        # Ensure regions begin at right STRIDE
        for region_i, region_data in self.g.regions.items():
            number_comparisons += 1
            expected_start = 1 + (region_data["ichr"] * STRIDE)
            self.assertEqual(expected_start, region_data["pos_start"])

            # -1 as the region includes the start element
            expected_end = (expected_start - 1) + LENGTH
            self.assertEqual(expected_end, region_data["pos_end"])
        self.assertEqual(self.EXPECTED_REGION_COUNT, number_comparisons)

###############################################################################
# \/ NOTE Tests below depend only on GroupOne                              \/ #
###############################################################################

    def test_content_counters(self):
        number_comparisons = 0
        for region_i, region_data in self.g.regions.items():
            number_comparisons += 1
            self.assertEqual(self.EXPECTED_REGIONS[region_data["chr"]][region_data["ichr"]], self.g.counter_matrix[self.g._get_group_id("GroupOne"), self.g._get_track_id("N"), region_i])
        self.assertEqual(self.EXPECTED_REGION_COUNT, number_comparisons)

    def test_number_buckets(self):
        self.assertEqual(len(self.EXPECTED_BUCKETS), len(self.g.group_buckets["GroupOne"]["N"]))
        total_in_bucket = 0
        total_in_bucket_g = 0

        for b in self.EXPECTED_BUCKETS:
            total_in_bucket += len(self.EXPECTED_BUCKETS[b])
        for b in self.g.group_buckets["GroupOne"]["N"]:
            total_in_bucket_g += len(self.g.group_buckets["GroupOne"]["N"][b])

        self.assertEqual(total_in_bucket, total_in_bucket_g)

    def test_content_buckets(self):
        number_comparisons = 0
        for bucket, content in self.g.group_buckets["GroupOne"]["N"].items():
            number_comparisons += 1
            self.assertEqual(sorted(self.EXPECTED_BUCKETS[bucket]), sorted(content))
        self.assertEqual(len(self.EXPECTED_BUCKETS), number_comparisons)

###############################################################################
# /\ NOTE Tests above depend only on GroupOne                              /\ #
###############################################################################

    def test_sort_min(self):
        number_comparisons = 0
        regions_min = self.g.query("min").candidates

        # Merge sort required for stability
        EXPECTED_MIN_REGION_ORDER = np.argsort(self.EXPECTED_REGION_TOTALS, kind="mergesort")

        for i, region in enumerate(regions_min):
            self.assertEquals(EXPECTED_MIN_REGION_ORDER[i], region['id'])
            number_comparisons += 1
        self.assertEqual(self.EXPECTED_REGION_COUNT, number_comparisons)


class TestGoldilocksRegression_NCounter_FASTA(unittest.TestCase):

#  GroupOne
#  CHR 1
#
#       1 |======***============***============***=========| 48
#          AAACCCNNNTTTAAACCCGGGNNNAAACCCGGGTTTNNNCCCGGGTTT
#          |0| |1| |1| |0| |0| |2| |0| |0| |0| |3| |0| |0|
#           |0| |2| |0| |0| |0| |3| |0| |0| |0| |2| |0| |0|
#            |0| |3| |0| |0| |0| |2| |0| |0| |1| |1| |0|
#             |0| |2| |0| |0| |1| |1| |0| |0| |2| |0| |0|
#

#  GroupOne
#  CHR 2
#
#       1 |====*=**==*=====*====*===**=====*=**===*==*===**| 48
#          GCGTNANNGANGGCTANTCTANAGCNNTTTCTNTNNGCANCANTTGNN
#          |0| |2| |1| |0| |1| |1| |2| |0| |2| |0| |1| |1|
#           |0| |2| |1| |0| |0| |1| |2| |0| |2| |1| |1| |2|
#            |1| |2| |1| |1| |0| |0| |1| |1| |2| |1| |1|
#             |1| |1| |0| |1| |1| |1| |0| |1| |1| |1| |0|
#
    @classmethod
    def setUpClass(cls):
        cls.g = Goldilocks(NucleotideCounterStrategy(["N"]), DATA_FAI, length=LENGTH, stride=STRIDE, is_faidx=True)
        cls.EXPECTED_REGIONS = {
            1: {
                0: 0,
                1: 0,
                2: 0,
                3: 0,
                4: 1,
                5: 2,
                6: 3,
                7: 2,
                8: 1,
                9: 0,
                10: 0,
                11: 0,
                12: 0,
                13: 0,
                14: 0,
                15: 0,
                16: 0,
                17: 0,
                18: 0,
                19: 1,
                20: 2,
                21: 3,
                22: 2,
                23: 1,
                24: 0,
                25: 0,
                26: 0,
                27: 0,
                28: 0,
                29: 0,
                30: 0,
                31: 0,
                32: 0,
                33: 0,
                34: 1,
                35: 2,
                36: 3,
                37: 2,
                38: 1,
                39: 0,
                40: 0,
                41: 0,
                42: 0,
                43: 0,
                44: 0,
                45: 0,
            },
            2: {
                0: 0,
                1: 0,
                2: 1,
                3: 1,
                4: 2,
                5: 2,
                6: 2,
                7: 1,
                8: 1,
                9: 1,
                10: 1,
                11: 0,
                12: 0,
                13: 0,
                14: 1,
                15: 1,
                16: 1,
                17: 0,
                18: 0,
                19: 1,
                20: 1,
                21: 1,
                22: 0,
                23: 1,
                24: 2,
                25: 2,
                26: 1,
                27: 0,
                28: 0,
                29: 0,
                30: 1,
                31: 1,
                32: 2,
                33: 2,
                34: 2,
                35: 1,
                36: 0,
                37: 1,
                38: 1,
                39: 1,
                40: 1,
                41: 1,
                42: 1,
                43: 0,
                44: 1,
                45: 2,
            }
        }

        cls.EXPECTED_REGION_TOTALS = []
        for j in cls.EXPECTED_REGIONS:
            for i in cls.EXPECTED_REGIONS[j]:
                cls.EXPECTED_REGION_TOTALS.append(cls.EXPECTED_REGIONS[j][i]*2)

        cls.EXPECTED_BUCKETS = {
            0: [0, 1, 2, 3, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 24, 25, 26,
                    27, 28, 29, 30, 31, 32, 33, 39, 40, 41, 42, 43, 44, 45,
                        46, 47, 57, 58, 59, 63, 64, 68, 73, 74, 75, 82, 89],
            1: [4, 8, 19, 23, 34, 38,
                        48, 49, 53, 54, 55, 56, 60, 61, 62, 65, 66, 67, 69,
                        72, 76, 77, 81, 83, 84, 85, 86, 87, 88, 90],
            2: [5, 7, 20, 22, 35, 37,
                        50, 51, 52, 70, 71, 78, 79, 80, 91],
            3: [6, 21, 36],
        }
        cls.EXPECTED_REGION_COUNT = len(cls.EXPECTED_REGIONS[1]) + len(cls.EXPECTED_REGIONS[2])

    def test_number_regions(self):
        self.assertEqual(self.EXPECTED_REGION_COUNT, len(self.g.regions.items()))

    # TODO Technically these check the metadata is correct and not the *actual* regions
    def test_region_lengths(self):
        number_comparisons = 0
        # Ensure ALL meet LENGTH
        for region_i, region_data in self.g.regions.items():
            number_comparisons += 1
            self.assertEqual(LENGTH, len(range(region_data["pos_start"], region_data["pos_end"])) + 1)
        self.assertEqual(self.EXPECTED_REGION_COUNT, number_comparisons)

    # TODO Technically these check the metadata is correct and not the *actual* regions
    def test_region_stride(self):
        number_comparisons = 0
        # Ensure regions begin at right STRIDE
        for region_i, region_data in self.g.regions.items():
            number_comparisons += 1
            expected_start = 1 + (region_data["ichr"] * STRIDE)
            self.assertEqual(expected_start, region_data["pos_start"])

            # -1 as the region includes the start element
            expected_end = (expected_start - 1) + LENGTH
            self.assertEqual(expected_end, region_data["pos_end"])
        self.assertEqual(self.EXPECTED_REGION_COUNT, number_comparisons)

###############################################################################
# \/ NOTE Tests below depend only on GroupOne                              \/ #
###############################################################################

    def test_content_counters(self):
        number_comparisons = 0
        for region_i, region_data in self.g.regions.items():
            number_comparisons += 1
            self.assertEqual(self.EXPECTED_REGIONS[region_data["chr"]][region_data["ichr"]], self.g.counter_matrix[self.g._get_group_id("GroupOne"), self.g._get_track_id("N"), region_i])
        self.assertEqual(self.EXPECTED_REGION_COUNT, number_comparisons)

    def test_number_buckets(self):
        self.assertEqual(len(self.EXPECTED_BUCKETS), len(self.g.group_buckets["GroupOne"]["N"]))
        total_in_bucket = 0
        total_in_bucket_g = 0

        for b in self.EXPECTED_BUCKETS:
            total_in_bucket += len(self.EXPECTED_BUCKETS[b])
        for b in self.g.group_buckets["GroupOne"]["N"]:
            total_in_bucket_g += len(self.g.group_buckets["GroupOne"]["N"][b])

        self.assertEqual(total_in_bucket, total_in_bucket_g)

    def test_content_buckets(self):
        number_comparisons = 0
        for bucket, content in self.g.group_buckets["GroupOne"]["N"].items():
            number_comparisons += 1
            self.assertEqual(sorted(self.EXPECTED_BUCKETS[bucket]), sorted(content))
        self.assertEqual(len(self.EXPECTED_BUCKETS), number_comparisons)

###############################################################################
# /\ NOTE Tests above depend only on GroupOne                              /\ #
###############################################################################

    def test_sort_min(self):
        number_comparisons = 0
        regions_min = self.g.query("min").candidates

        # Merge sort required for stability
        EXPECTED_MIN_REGION_ORDER = np.argsort(self.EXPECTED_REGION_TOTALS, kind="mergesort")

        for i, region in enumerate(regions_min):
            self.assertEquals(EXPECTED_MIN_REGION_ORDER[i], region['id'])
            number_comparisons += 1
        self.assertEqual(self.EXPECTED_REGION_COUNT, number_comparisons)

if __name__ == '__main__':
    unittest.main()
