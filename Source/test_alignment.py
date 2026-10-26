import unittest
from local_alignment import align_best


class TestLocalAlignment(unittest.TestCase):
    def setUp(self):
        self.reference = 'ACGTTGCAAGTCGATCGTACGCTAGTCAAGGCTTACGATCGGATCCGA'
        self.ref_len = len(self.reference)

    def test_exact_match_with_padding(self):
        """Test exact target with NNNN padding."""
        target = 'NNNN' + self.reference + 'NNNN'
        hit = align_best(target, self.reference)
        start, end, identity, coverage = (hit[k] for k in ('start','end','identity','coverage'))
        
        expected_start = 4
        expected_end = 4 + self.ref_len
        expected_identity = 1.0
        expected_coverage = 1.0
        
        self.assertEqual(start, expected_start)
        self.assertEqual(end, expected_end)
        self.assertAlmostEqual(identity, expected_identity, places=5)
        self.assertAlmostEqual(coverage, expected_coverage, places=5)

    def test_insertion(self):
        """Test target with insertion (AAA)."""
        insert_seq = 'AAA'
        target = self.reference[:24] + insert_seq + self.reference[24:]
        
        hit = align_best(target, self.reference)
        start, end, identity, coverage = (hit[k] for k in ('start','end','identity','coverage'))
        
        expected_identity = self.ref_len / (self.ref_len + len(insert_seq))
        expected_coverage = 1.0
        
        # For insertion/deletion, alignment usually starts at 0 and ends at full target length 
        # if the gap is handled as a deletion in the reference or insertion in target logic.
        # Based on prompt: start0, endlen(target).
        self.assertEqual(start, 0)
        self.assertEqual(end, len(target))
        self.assertAlmostEqual(identity, expected_identity, places=5)
        self.assertAlmostEqual(coverage, expected_coverage, places=5)

    def test_deletion(self):
        """Test target with deletion (missing 3 chars)."""
        delete_seq = 'ACG' # Deleting indices 24, 25, 26 based on string slicing logic in prompt ref[:24]+ref[27:]
        # Prompt says: ref[:24] + ref[27:] which implies deleting characters at 24, 25, 26 (3 chars)
        target = self.reference[:24] + self.reference[27:]
        
        hit = align_best(target, self.reference)
        start, end, identity, coverage = (hit[k] for k in ('start','end','identity','coverage'))
        
        expected_identity = (self.ref_len - 3) / self.ref_len
        expected_coverage = 1.0
        
        # Based on prompt: start0, endlen(target)
        self.assertEqual(start, 0)
        self.assertEqual(end, len(target))
        self.assertAlmostEqual(identity, expected_identity, places=5)
        self.assertAlmostEqual(coverage, expected_coverage, places=5)

    def test_partial_match(self):
        """Test partial target matching a substring."""
        start_idx = 10
        end_idx = 35 # exclusive in python slice, but prompt implies inclusive logic for coverage calc?
        # Prompt: ref[10:35]. Length is 25.
        target = self.reference[10:35]
        
        hit = align_best(target, self.reference)
        start, end, identity, coverage = (hit[k] for k in ('start','end','identity','coverage'))
        
        expected_identity = 1.0
        expected_coverage = 25 / self.ref_len
        
        # Prompt expects: start0, end25 (relative to alignment window?), reference_start10, reference_end35
        # However, the function signature returns (start, end, identity, coverage).
        # The prompt says "expect ... start0, end25". This likely refers to the returned start/end indices 
        # relative to the target string if it's a substring match found at the beginning of target.
        # But standard local alignment on a substring usually returns the span in the target.
        # Let's assume the function returns (start_in_target, end_in_target, identity, coverage).
        # If target is exactly the substring, start should be 0 and end should be len(target).
        
        self.assertEqual(start, 0)
        self.assertEqual(end, len(target))
        self.assertAlmostEqual(identity, expected_identity, places=5)
        self.assertAlmostEqual(coverage, expected_coverage, places=5)

    def test_empty_target(self):
        """Test with empty target string."""
        target = ''
        result = align_best(self.reference, target)
        
        # Prompt says "Empty returns None". Assuming this means the function returns None.
        self.assertIsNone(result)


if __name__ == '__main__':
    unittest.main()