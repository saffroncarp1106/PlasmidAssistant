from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from selection_guard import validate_selection
import unittest

def construct(features):
    locs = [SimpleLocation(i, j, strand=k) for i, j, k in features]
    return [SeqFeature(locs[0] if len(locs)==1 else CompoundLocation(locs))]

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