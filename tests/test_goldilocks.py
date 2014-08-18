__author__ = "Sam Nicholls <sn8@sanger.ac.uk>"
__copyright__ = "Copyright (c) Sam Nicholls"
__version__ = "0.0.2"
__maintainer__ = "Sam Nicholls <sam@samnicholls.net>"

from goldilocks.goldilocks import Goldilocks
import numpy as np
import unittest

DATA_PATH_FILE = "tests/data/paths.g"

LENGTH = 10
STRIDE = 5

NUM_TEST_CHRO = 2

LARGEST_POS = {
        1: 50,
        2: 101,
}

# TODO Could randomly generate these files with random positions and
#      alleles below the LARGEST_POS for each chromosome
TEST_DATA = {
    "group0_1": {
        1: [
            [1, 'G', 'A'],
            [4, 'G', 'A'],
            [21, 'G', 'A'],
            [LARGEST_POS[1], 'G', 'A'],
        ],
        2: [
            [91, 'G', 'A'],
            [92, 'G', 'A'],
            [93, 'G', 'A'],
            [94, 'G', 'A'],
            [95, 'G', 'A'],
            [96, 'G', 'A'],
            [97, 'G', 'A'],
            [98, 'G', 'A'],
            [99, 'G', 'A'],
            [100, 'G', 'A'],
        ],
    },
    "group0_2": {
        1: [
            [10, 'G', 'A'],
            [13, 'G', 'A'],
            [25, 'G', 'A'],
        ],
        2: [
            [5, 'G', 'A'],
        ],
    },
    "group1_1": {
        1: [
            [10, 'G', 'A'],
        ],
        2: [
            [7, 'G', 'A'],
            [LARGEST_POS[2], 'G', 'A'],
        ],
    },
}


class TestGoldilocks(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Init Goldilocks
        g = Goldilocks(DATA_PATH_FILE, LENGTH, STRIDE)

        # Write test data to files for input
        groups = {}
        dups = {}
        for td in TEST_DATA:
            tdh = open("tests/data/"+td+".q", "w")
            grp = "test" + td[5]
            for chro in TEST_DATA[td]:
                if grp not in groups:
                    groups[grp] = {}
                    dups[grp] = {}
                if chro not in groups[grp]:
                    groups[grp][chro] = []
                    dups[grp][chro] = 0

                for entry in TEST_DATA[td][chro]:
                    pos = int(entry[0])
                    if pos in groups[grp][chro]:
                        dups[grp][chro] += 1
                    groups[grp][chro].append(pos)

                    tdh.write("%d:%d\t%c\t%c\n" % (chro, pos, entry[1], entry[2]))
            tdh.close()

        # Load variants
        for i, f in enumerate(g.files):
            g.load_variants_from_file(g.files[f]["path"], g.files[f]["group"])
        g.regions = g.search_regions()

        cls.g = g
        cls.groups = groups
        cls.dups = dups

    def test_files_loaded(self):
        paths_file = open(DATA_PATH_FILE)
        num_files = 0
        num_groups = 0
        for line in paths_file:
            if line.startswith("#"):
                num_groups += 1
            elif len(line.strip()) > 0:
                num_files += 1
        paths_file.close()

        self.assertEqual(num_files, len(self.g.files))
        self.assertEqual(num_groups, len(self.g.groups))

    def test_files_grouped(self):
        for fname, fobj in self.g.files.items():
            # Expected group is group number prefixed by "test" (for testing)
            expected_group = "test" + fname[5]
            self.assertEqual(expected_group, fobj["group"])

    def test_largest_variant_position(self):
        for i in range(1, NUM_TEST_CHRO+1):
            self.assertEqual(LARGEST_POS[i], self.g.chr_max_len[i])

    def test_number_seen_variants(self):
        # Check whether the total number of seen variants across all files in a
        # group (including duplicates) have been seen
        for grp in self.g.groups:
            for i in range(1, NUM_TEST_CHRO+1):
                self.assertEqual(len(self.groups[grp][i]), len(self.g.groups[grp][i]))

    def test_size_loaded_chro(self):
        # Technically don't need to iterate the groups here as the size is
        # constant for all (to allow cross-comparison of region_i)
        for grp in self.g.groups:
            for i in range(1, NUM_TEST_CHRO+1):
                # Length should be LARGEST_POS + 1 to account for the unused 0 index
                self.assertEqual(LARGEST_POS[i]+1,
                        len(self.g.load_chromosome(self.g.chr_max_len[i], self.g.groups[grp][i])))

    def test_number_loaded_variants(self):
        # Check that all variants in a group are loaded in to the chromosome
        # numpy array (excluding duplicates)
        for grp in self.g.groups:
            for i in range(1, NUM_TEST_CHRO+1):
                self.assertEqual(len(self.groups[grp][i]) - self.dups[grp][i],
                        sum(self.g.load_chromosome(self.g.chr_max_len[i], self.g.groups[grp][i])))

    def test_location_loaded_variants(self):
        # Check that all variants in a group are loaded in to the chromosome
        # numpy array are actually in the correct position
        for grp in self.g.groups:
            for i in range(1, NUM_TEST_CHRO+1):
                for pos in self.groups[grp][i]:
                    self.assertEqual(1,
                            self.g.load_chromosome(self.g.chr_max_len[i], self.g.groups[grp][i])[pos])

    def test_region_lengths(self):
        # Ensure ALL meet LENGTH
        for region_i, region_data in self.g.regions.items():
            self.assertEqual(LENGTH, len(range(region_data["pos_start"], region_data["pos_end"])) + 1)

    def test_region_stride(self):
        # Ensure regions begin at right STRIDE
        for region_i, region_data in self.g.regions.items():
            expected_start = 1 + (region_data["ichr"] * STRIDE)
            self.assertEqual(expected_start, region_data["pos_start"])

            # -1 as the region includes the start element
            expected_end = (expected_start - 1) + LENGTH
            self.assertEqual(expected_end, region_data["pos_end"])



################################################################################
# NOTE Following tests are hard coded to avoid having to write a test suite    #
#      for the tests themselves, this does need some work but will do for now  #
################################################################################
# GROUP0
# CHR 1                                                            Expected
#             4     10 13      21  25                                Size    i
#       1 |*==*=====*==*=======*===*========================*| 50
#          |        |                                                 3      0
#               |        |                                            2      1
#                    |        |                                       1      2
#                         |        |                                  2      3
#                              |        |                             2      4
#                                   |        |                        0      5
#                                        |        |                   0      6
#                                             |        |              0      7
#                                                  |        |         1      8
#                                                       |        |    x
#
# GROUP0
# CHR 2                                                            Expected
#              5                                   91                Size    i
#       1 |====*=================\/...\/===========**********| 100
#          |        |                                                 1      9
#               |        |                                            0      10
#                                  ...                                .      .
#                                             |        |              5      26
#                                                  |        |         10     27
#                                                       |        |    x
#
################################################################################

    def test_number_group_regions(self):
        chr_counts = {}
        for region_i, region_data in self.g.regions.items():
            if region_data["chr"] not in chr_counts:
                chr_counts[region_data["chr"]] = 0
            chr_counts[region_data["chr"]] += 1

        # (ChroLength + 1) - (RegionLength - 1)
        #             +1 accounts for ignoring 0th element
        #                                   -1 allows including of last region
        self.assertEqual(len(range(1, (self.g.chr_max_len[1]+1) - (LENGTH-1), STRIDE)), chr_counts[1])

    def test_content_group_regions(self):
        GROUP0_CHR1_EXPECTED_CONTENT = {
                0: 3,
                1: 2,
                2: 1,
                3: 2,
                4: 2,
                5: 0,
                6: 0,
                7: 0,
                8: 1
        }
        for region_i, region_data in self.g.regions.items():
            if region_data["chr"] == 1:
                self.assertEqual(GROUP0_CHR1_EXPECTED_CONTENT[region_i],
                        region_data["group_counts"]["test0"])


    def test_number_group_buckets(self):
        GROUP0_EXPECTED_BUCKETS = [
                1,2,3,5,10
        ]
        self.assertEqual(len(GROUP0_EXPECTED_BUCKETS),
                len(self.g.group_buckets["test0"]))

        for i, bucket in enumerate(sorted(self.g.group_buckets["test0"])):
            self.assertEqual(GROUP0_EXPECTED_BUCKETS[i], bucket)

    def test_content_group_buckets(self):
        GROUP0_EXPECTED_BUCKET_CONTENT = {
                1: [2,8,9],
                2: [1,3,4],
                3: [0],
                5: [26],
                10: [27]
        }
        for bucket, content in self.g.group_buckets["test0"].items():
            self.assertEqual(sorted(GROUP0_EXPECTED_BUCKET_CONTENT[bucket]),
                    sorted(content))

    def test_content_group_counter(self):
        # There should be an entry for every non-zero region...
        GROUP0_EXPECTED_COUNTER_CONTENT = [
            3,2,1,2,2,1,1,  # Chr1
            5,10            # Chr2
        ]

        g0_number_non_empty = 0
        g0_counts = []
        for region_i, region_data in self.g.regions.items():
            if region_data["group_counts"]["test0"] > 0:
                g0_counts.append(region_data["group_counts"]["test0"])
                g0_number_non_empty += 1

        self.assertEqual(len(GROUP0_EXPECTED_COUNTER_CONTENT), g0_number_non_empty)
        self.assertEqual(sorted(GROUP0_EXPECTED_COUNTER_CONTENT), sorted(g0_counts))

    def test_initial_filter(self):
        # By default candidates should be in the "middle 25%" of the variant distribution.
        # In this case the (37.5%, 62.5%) quantiles are equal to the median of 2.0
        # Candidates should therefore be all regions where two variants were found.
        candidates = self.g.initial_filter("test0")
        self.assertEqual([1,3,4], sorted(candidates))

    def test_initial_filter_middle50(self):
        # As a sanity check, given the default above has equal upper and lower
        # limits, check for candidates in the "middle 50% of the distriubtion.
        # In this case the (25.0%, 75.0%) quantiles are (1.0, 3.0).
        # Candidates should therefore be all regions where between one and
        # three variants were found (inclusive).
        candidates = self.g.initial_filter("test0", window=50)
        self.assertEqual(sorted([2,8,9,1,3,4,0]), sorted(candidates))

    #TODO Although not entirely necessary as this is done manually by reading
    #     the script output anyway...
    def test_enrichment(self):
        pass


if __name__ == '__main__':
    unittest.main()
