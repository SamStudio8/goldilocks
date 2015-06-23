import unittest

import numpy as np

from goldilocks.goldilocks import Goldilocks
from goldilocks.strategies import StrategyValue
from goldilocks.strategies import BaseStrategy, NucleotideCounterStrategy

################################################################################
# NOTE Tests following do not test the correctness of regions located by
#      Goldilocks but are instead designed to ensure the program correctly
#      handles inputs and outputs as expected.
################################################################################

sequence_data = {
        "my_sample": {
            "one": "CATCANCAT",
            2: "NANANANANA",
            "X": "GATTACAGATTACAN",
            "three": "..A",
        },
        "my_other_sample": {
            "one": "TATANTATA",
            2: "GANGANGAN",
            "X": "GATTACAGATTACAN",
            "three": ".N.",
        }
}

OPS = ["median", "max", "min", "mean"]

class TestGoldilocks(unittest.TestCase):

    def setUp(self):
        self.g = Goldilocks(NucleotideCounterStrategy(["A","C","G","T","N"]), sequence_data, length=3, stride=1)
        self.TOTAL_REGIONS = 29

    def __test_simple_exclusions(self, EXCLUSIONS, limit=0):

        FILTER_TO_PROPERTY = {
            "start_lte": ("pos_start", "lt"),
            "start_gte": ("pos_start", "gt"),
            "end_lte": ("pos_end", "lt"),
            "end_gte": ("pos_end", "gt"),
            "chr": ("chr", "nin")
        }

        for exclusion_name, exclusion in EXCLUSIONS.items():
            for op in OPS:
                if limit > 0:
                    candidates = self.g._filter(op, exclusions={
                        exclusion["filter"]: exclusion["value"]
                    }, limit=limit).candidates
                else:
                    candidates = self.g._filter(op, exclusions={
                        exclusion["filter"]: exclusion["value"]
                    }).candidates

                for c in candidates:
                    cproperty = FILTER_TO_PROPERTY[exclusion["filter"]][0]
                    test_type = FILTER_TO_PROPERTY[exclusion["filter"]][1]

                    if test_type == "lt":
                        self.assertTrue(c[cproperty] > exclusion["value"])
                    elif test_type == "gt":
                        self.assertTrue(c[cproperty] < exclusion["value"])
                    elif test_type == "nin":
                        self.assertNotIn(c[cproperty], exclusion["value"])
                    else:
                        self.fail("Incorrect test_type")

                if len(candidates) == 0 and ("expect_none" not in exclusion):
                    self.fail("No candidates returned but at least one expected...")

                if limit:
                    if limit > self.TOTAL_REGIONS:
                        # Don't test if limit is larger than number of regions
                        pass
                    elif "expect_none" not in exclusion:
                        self.assertEqual(limit, len(candidates))
                    else:
                        self.assertEqual(0, len(candidates))

    def test_missing_length(self):
        self.assertRaises(TypeError, Goldilocks, NucleotideCounterStrategy([]), sequence_data, stride=1)

    def test_missing_stride(self):
        self.assertRaises(TypeError, Goldilocks, NucleotideCounterStrategy([]), sequence_data, length=1)

    def test_invalid_stride(self):
        self.assertRaises(ValueError, Goldilocks, NucleotideCounterStrategy([]), sequence_data, length=1, stride=0)
        self.assertRaises(ValueError, Goldilocks, NucleotideCounterStrategy([]), sequence_data, length=1, stride=-1)
        self.assertRaises(ValueError, Goldilocks, NucleotideCounterStrategy([]), sequence_data, length=1, stride=-1000)

    def test_invalid_length(self):
        self.assertRaises(ValueError, Goldilocks, NucleotideCounterStrategy([]), sequence_data, length=0, stride=1)
        self.assertRaises(ValueError, Goldilocks, NucleotideCounterStrategy([]), sequence_data, length=-1, stride=1)
        self.assertRaises(ValueError, Goldilocks, NucleotideCounterStrategy([]), sequence_data, length=-1000, stride=1)

    def test_invalid_filter_distance(self):
        for op in OPS:
            self.assertRaises(ValueError, self.g._filter, op, actual_distance=1, percentile_distance=1)

    def test_invalid_sort_operation(self):
        for op in OPS:
            self.assertRaises(TypeError, self.g._filter, "hoot")

    def test_unimplemented_strategy(self):
        self.assertRaises(NotImplementedError, Goldilocks, BaseStrategy(), sequence_data, length=1, stride=1)

    def test_exclude_chr(self):
        EXCLUSIONS = {
            "simple_chr_str": {
                "filter": "chr",
                "value": ["one"],
            },
            "simple_chr_int": {
                "filter": "chr",
                "value": [2],
            },
            "simple_all_chr": {
                "filter": "chr",
                "value": ["one", "X", 2, "three"],
                "expect_none": True
            },
        }

        self.__test_simple_exclusions(EXCLUSIONS)
        self.__test_simple_exclusions(EXCLUSIONS, limit=1)
        self.__test_simple_exclusions(EXCLUSIONS, limit=5)
        self.__test_simple_exclusions(EXCLUSIONS, limit=100)

    def test_exclude_start_gte(self):
        EXCLUSIONS = {
            "simple_start_gte": {
                "filter": "start_gte",
                "value": 5,
            },
            "prevent_start_gte": {
                "filter": "start_gte",
                "value": 100,
            },
            "all_start_gte": {
                "filter": "start_gte",
                "value": 1,
                "expect_none": True
            },
        }

        self.__test_simple_exclusions(EXCLUSIONS)
        self.__test_simple_exclusions(EXCLUSIONS, limit=1)
        self.__test_simple_exclusions(EXCLUSIONS, limit=5)
        self.__test_simple_exclusions(EXCLUSIONS, limit=100)


    def test_exclude_start_lte(self):
        EXCLUSIONS = {
            "simple_start_lte": {
                "filter": "start_lte",
                "value": 5,
            },
            "prevent_start_lte": {
                "filter": "start_lte",
                "value": 0,
            },
            "all_start_lte": {
                "filter": "start_lte",
                "value": 100,
                "expect_none": True
            },
        }

        self.__test_simple_exclusions(EXCLUSIONS)
        self.__test_simple_exclusions(EXCLUSIONS, limit=1)
        self.__test_simple_exclusions(EXCLUSIONS, limit=5)
        self.__test_simple_exclusions(EXCLUSIONS, limit=100)

    def test_exclude_end_gte(self):
        EXCLUSIONS = {
            "simple_end_gte": {
                "filter": "end_gte",
                "value": 5,
            },
            "prevent_end_gte": {
                "filter": "end_gte",
                "value": 100,
            },
            "all_end_gte": {
                "filter": "end_gte",
                "value": 1,
                "expect_none": True
            },
        }

        self.__test_simple_exclusions(EXCLUSIONS)
        self.__test_simple_exclusions(EXCLUSIONS, limit=1)
        self.__test_simple_exclusions(EXCLUSIONS, limit=5)
        self.__test_simple_exclusions(EXCLUSIONS, limit=100)

    def test_exclude_end_lte(self):
        EXCLUSIONS = {
            "simple_end_lte": {
                "filter": "end_lte",
                "value": 5,
            },
            "prevent_end_lte": {
                "filter": "end_lte",
                "value": 0,
            },
            "all_end_lte": {
                "filter": "end_lte",
                "value": 100,
                "expect_none": True
            },
        }

        self.__test_simple_exclusions(EXCLUSIONS)
        self.__test_simple_exclusions(EXCLUSIONS, limit=1)
        self.__test_simple_exclusions(EXCLUSIONS, limit=5)
        self.__test_simple_exclusions(EXCLUSIONS, limit=100)

    def test_exclude_and(self):
        for op in OPS:
            candidates = self.g._filter(op, exclusions={
                                                        "start_gte": 5,
                                                        "end_lte": 9,
                                                        }, use_and=True).candidates
            for c in candidates:
                self.assertFalse(c["pos_start"] >= 5 and c["pos_end"] <= 9)

            self.assertTrue(len(candidates) < self.TOTAL_REGIONS)

    def test_exclude_and_with_chr(self):
        for op in OPS:
            candidates = self.g._filter(op, exclusions={
                                                        "start_gte": 5,
                                                        "end_lte": 9,
                                                        "chr": ["X"],
                                                        }, use_and=True).candidates
            non_x_count = 0
            for c in candidates:
                if c["chr"] == "X":
                    self.assertFalse(c["pos_start"] >= 5 and c["pos_end"] <= 9)
                else:
                    if (c["pos_start"] >= 5 and c["pos_end"] <= 9):
                        non_x_count += 1
            self.assertTrue(non_x_count > 0)

            self.assertTrue(len(candidates) < self.TOTAL_REGIONS)

    def test_exclude_chr_specific_chr(self):
        pass

    def test_exclude_chr_specific_start(self):
        pass

    def test_exclude_chr_specific_end(self):
        pass

    def test_exclude_chr_specific_and(self):
        pass

    def test_limit(self):
        for op in OPS:
            candidates = self.g.query(op, limit=1).candidates
            self.assertTrue(len(candidates) == 1)

            candidates = self.g.query(op, limit=10).candidates
            self.assertTrue(len(candidates) == 10)

            candidates = self.g.query(op, limit=100).candidates
            self.assertTrue(len(candidates) == self.TOTAL_REGIONS)

    def test_distance_upper(self):
        pass

    def test_distance_lower(self):
        pass

    def test_distance_around(self):
        pass

class TestStrategyValue(unittest.TestCase):

    def test_default_k_0(self):
        number_comparisons = 0
        values = [0, 0.1, 1.0, -1.0, 10.0, 100.0]
        for value in values:
            sv = StrategyValue(value)
            self.assertEquals(value, sv)
            self.assertEquals(0.0, sv.k)
            number_comparisons += 1
        self.assertEqual(len(values), number_comparisons)

    def test_k(self):
        number_comparisons = 0
        values = [0.1, 1.0, -1.0, 10.0, 100.0, 0, 0.0]
        weights = [0, 1, 10, 100]
        for value in values:
            for k in weights:
                sv = StrategyValue(value, k)
                self.assertEquals(value, sv)
                self.assertEquals(k, sv.k)
                number_comparisons += 1
        self.assertEqual(len(values) * len(weights), number_comparisons)

    def test_invalid_k(self):
        number_comparisons = 0
        values = [0.1, 1.0, -1.0, 10.0, 100.0, 0, 0.0]
        weights = [-1, -1.0, -0.1, -100, -10.0, 0.1, 0.5, 0.999]
        for value in values:
            for k in weights:
                self.assertRaises(ValueError, StrategyValue, value, k)
                number_comparisons += 1
        self.assertEqual(len(values) * len(weights), number_comparisons)

    def test_weighted_add(self):
        cascade_3 = StrategyValue(82, 180)
        cascade_2 = StrategyValue(70.25, 120)
        cascade_1 = StrategyValue(72.5, 60)
        self.assertEqual(76.5, 0 + cascade_3 + cascade_2 + cascade_1)
        self.assertEqual(76.5, cascade_3 + cascade_2 + cascade_1 + 0)
        self.assertEqual(76.5, 1 + cascade_3 + cascade_2 + cascade_1 + (-1))
        self.assertEqual(76.5, StrategyValue(1) + cascade_3 + cascade_2 + cascade_1 + (-1))
        self.assertEqual(76.5, StrategyValue(1) + cascade_3 + cascade_2 + cascade_1 + StrategyValue(-1))
        self.assertEqual(76.5, StrategyValue(0) + cascade_3 + cascade_2 + cascade_1)

        values_a = [-10, -1, 0, 1, 3, 5, 7, 9, 15, 50, 100, 1000]
        weights_a = [1, 3, 5, 10, 25, 100]
        values_b = [-50, -5, 0, 0.1, 11, 22, 33, 44, 55, 1000000]
        weights_b = [1, 2, 9, 15, 50, 1000]

        number_comparisons = 0
        for a in values_a:
            for b in values_b:
                for w_a in weights_a:
                    for w_b in weights_b:
                        sv = StrategyValue(a, w_a) + StrategyValue(b, w_b)
                        total_a = a * w_a
                        total_b = b * w_b
                        expected = (total_a + total_b) / float(w_a + w_b)
                        self.assertEqual(expected, sv)
                        number_comparisons += 1
        self.assertEqual(len(values_a) * len(values_b) * len(weights_a) * len(weights_b), number_comparisons)

    def test_weight(self):
        self.assertEqual(-49.9, 1 + StrategyValue(-50, 10))
        self.assertEqual(50, StrategyValue(0, 100) + StrategyValue(100, 100))

        self.assertEqual(2, StrategyValue(1, 0) + StrategyValue(1, 0))
        self.assertEqual(2, 1 + StrategyValue(1, 0))
        self.assertEqual(2, StrategyValue(1, 0) + 1)

        self.assertEqual(1, StrategyValue(0.5, 0) + StrategyValue(0.5, 0))
        self.assertEqual(1, 0.5 + StrategyValue(0.5, 0))
        self.assertEqual(1, StrategyValue(0.5, 0) + 0.5)

        self.assertEqual(2, 1 + StrategyValue(1, 1))
        self.assertEqual(2, StrategyValue(1, 1) + 1)
        self.assertEqual(2, StrategyValue(1, 0) + StrategyValue(1, 1))

        self.assertEqual(1, StrategyValue(1, 1) + StrategyValue(1, 1))
        self.assertEqual(2, StrategyValue(2, 1) + StrategyValue(2, 1))

        self.assertEqual(10.0, StrategyValue(0) + StrategyValue(10, 100))
        self.assertEqual(10.0, 0 + StrategyValue(10, 100))
        self.assertEqual(10.01, 1 + StrategyValue(10, 100))

    def test_weighted_radd(self):
        values_a = [-10, -1, 0, 1, 3, 5, 7, 9, 15, 50, 100, 1000]
        values_b = [-50, -5, 0, 0.1, 11, 22, 33, 44, 55, 1000000]
        weights_b = [1, 2, 9, 15, 50, 1000]

        number_comparisons = 0
        for a in values_a:
            for b in values_b:
                for w_b in weights_b:
                    sv = a + StrategyValue(b, w_b)
                    total_a = a
                    total_b = b * w_b
                    expected = (total_a + total_b) / float(0 + w_b)
                    self.assertEqual(expected, sv)
                    number_comparisons += 1
        self.assertEqual(len(values_a) * len(values_b) * len(weights_b), number_comparisons)

    def test_unweighted_add(self):
        number_comparisons = 0
        values_a = [-10, -1, 0, 1, 3, 5, 7, 9, 15, 50, 100, 1000]
        values_b = [-50, -5, 0, 0.1, 11, 22, 33, 44, 55, 1000000]
        for a in values_a:
            for b in values_b:
                sv = StrategyValue(a, 0) + StrategyValue(b, 0)
                self.assertEqual(float(a) + b, float(sv))
                number_comparisons += 1
        self.assertEqual(len(values_a) * len(values_b), number_comparisons)

    def test_non_sv_add(self):
        pass

    def test_non_sv_radd(self):
        pass


if __name__ == '__main__':
    unittest.main()
