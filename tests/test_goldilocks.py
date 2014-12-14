import unittest

import numpy as np

from goldilocks.goldilocks import Goldilocks
from goldilocks.strategies import NucleotideCounterStrategy

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

    @classmethod
    def setUpClass(cls):
        cls.g = Goldilocks(NucleotideCounterStrategy(["A","C","G","T","N"]), sequence_data, length=3, stride=1)
        cls.TOTAL_REGIONS = 29

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
                    }, limit=limit)
                else:
                    candidates = self.g._filter(op, exclusions={
                        exclusion["filter"]: exclusion["value"]
                    })

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

    def test_invalid_stride(self):
        self.assertRaises(Exception, Goldilocks, NucleotideCounterStrategy([]), sequence_data, stride=0)
        self.assertRaises(Exception, Goldilocks, NucleotideCounterStrategy([]), sequence_data, stride=-1)
        self.assertRaises(Exception, Goldilocks, NucleotideCounterStrategy([]), sequence_data, stride=-1000)

    def test_invalid_length(self):
        self.assertRaises(Exception, Goldilocks, NucleotideCounterStrategy([]), sequence_data, length=0)
        self.assertRaises(Exception, Goldilocks, NucleotideCounterStrategy([]), sequence_data, length=-1)
        self.assertRaises(Exception, Goldilocks, NucleotideCounterStrategy([]), sequence_data, length=-1000)

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
                                                        }, use_and=True)
            for c in candidates:
                self.assertFalse(c["pos_start"] >= 5 and c["pos_end"] <= 9)

            self.assertTrue(len(candidates) < self.TOTAL_REGIONS)

    def test_exclude_and_with_chr(self):
        for op in OPS:
            candidates = self.g._filter(op, exclusions={
                                                        "start_gte": 5,
                                                        "end_lte": 9,
                                                        "chr": ["X"],
                                                        }, use_and=True)
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
            candidates = self.g._filter(op, limit=1)
            self.assertTrue(len(candidates) == 1)

            candidates = self.g._filter(op, limit=10)
            self.assertTrue(len(candidates) == 10)

            candidates = self.g._filter(op, limit=100)
            self.assertTrue(len(candidates) == self.TOTAL_REGIONS)

if __name__ == '__main__':
    unittest.main()
