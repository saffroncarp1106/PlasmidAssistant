from Bio import SeqFeature, SeqLocation
import unittest

def construct(features):
    locs = [SeqFeature(SeqLocation(i, j, k), i) for i, j, k in features]
    return locs

class TestSelectionGuard(unittest.TestCase):
    def test_same_strand_overlap(self):
        f1 = construct([(0, 100, 1)])
        f2 = construct([(5, 105, 1)])
        with self.assertRaises(ValueError):
            validate_selection(f1 + f2)

    def test_opposite_strands_allowed(self):
        f1 = construct([(0, 100, 1)])
        f2 = construct([(5, 105, -1)])
        self.assertIsNone(validate_selection(f1 + f2))

    def test_circular_disjoint(self):
        f1 = construct([(90, 100, 1), (0, 10, 1)])
        f2 = construct([(20, 80, 1)])
        self.assertIsNone(validate_selection(f1 + f2))

if __name__ == '__main__':
    unittest.main()