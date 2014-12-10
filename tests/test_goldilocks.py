import unittest

import numpy as np

from goldilocks.goldilocks import Goldilocks
from goldilocks.strategies import NCounterStrategy

################################################################################
# NOTE Tests following do not test the correctness of regions located by
#      Goldilocks but are instead designed to ensure the program correctly
#      handles inputs and outputs as expected.
################################################################################

DATA = {
    "GroupOne": {
        1: "AAACCCNNNTTTAAACCCGGGNNNAAACCCGGGTTTNNNCCCGGGTTT",
        2: "GCGTNANNGANGGCTANTCTANAGCNNTTTCTNTNNGCANCANTTGNN",
    }
}

class TestGoldilocks(unittest.TestCase):

    def test_invalid_stride(self):
        self.assertRaises(Exception, Goldilocks, NCounterStrategy(), DATA, stride=0)
        self.assertRaises(Exception, Goldilocks, NCounterStrategy(), DATA, stride=-1)
        self.assertRaises(Exception, Goldilocks, NCounterStrategy(), DATA, stride=-1000)

    def test_invalid_length(self):
        self.assertRaises(Exception, Goldilocks, NCounterStrategy(), DATA, length=0)
        self.assertRaises(Exception, Goldilocks, NCounterStrategy(), DATA, length=-1)
        self.assertRaises(Exception, Goldilocks, NCounterStrategy(), DATA, length=-1000)


if __name__ == '__main__':
    unittest.main()
