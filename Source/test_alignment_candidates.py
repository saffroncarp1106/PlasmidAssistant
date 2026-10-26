import random
import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from alignment_candidates import alignment_candidates
from Bio.SeqFeature import SeqFeature, SimpleLocation

# Fixture setup
rng = random.Random(123)
REFERENCE = ''.join(rng.choice('ACGT') for _ in range(120))
ENTRY_FIELDS = {
    'sequence': REFERENCE,
    'feature_id': 'fixture',
    'canonical_name': 'Synthetic',
    'feature_type': 'misc_feature',
    'source': 'synthetic',
    'version': '1'
}

class TestAlignmentCandidates(unittest.TestCase):

    def test_partial_match(self):
        record = SeqRecord(Seq(REFERENCE[20:110]), id='test_ref')
        entry = ENTRY_FIELDS.copy()
        fast_hits = []
        
        candidates = sorted(alignment_candidates(record, [entry], fast_hits), key=lambda h: (h['coverage'],h['identity']), reverse=True)
        
        self.assertTrue(any(c['match'] == 'partial' for c in candidates))
        self.assertGreaterEqual(candidates[0]['identity'], 0.75)
        self.assertEqual(set(candidates[0]['location']), set(range(90)))
        self.assertNotIn('High confidence', [c['confidence'] for c in candidates])

    def test_insertion_reverse(self):
        target = REFERENCE[:60] + 'AAA' + REFERENCE[60:]
        record = SeqRecord(Seq(target).reverse_complement(), id='test_ref_ins')
        entry = ENTRY_FIELDS.copy()
        
        # Create a fast hit to trigger suppression logic if needed, but here we test reverse complement detection
        fast_hits = []
        
        candidates = sorted(alignment_candidates(record, [entry], fast_hits), key=lambda h: (h['coverage'],h['identity']), reverse=True)
        
        self.assertTrue(any(c['location'].strand == -1 for c in candidates))
        self.assertGreaterEqual(candidates[0]['coverage'], 0.95)
        extracted_seq = SeqFeature(candidates[0]['location']).extract(record.seq).upper()
        self.assertEqual(extracted_seq, target.upper())

    def test_circular_topology(self):
        target = REFERENCE[40:] + 'NNNN' + REFERENCE[:40]
        record = SeqRecord(Seq(target), id='test_ref_circ', annotations={'topology': 'circular'})
        entry = ENTRY_FIELDS.copy()
        
        fast_hits = []
        candidates = sorted(alignment_candidates(record, [entry], fast_hits), key=lambda h: (h['coverage'],h['identity']), reverse=True)
        
        self.assertTrue(any(c['location'].strand == 1 for c in candidates))
        # Check that location parts cover the circular shift logic implicitly via identity/coverage
        self.assertGreaterEqual(candidates[0]['identity'], 0.95)
        extracted_seq = SeqFeature(candidates[0]['location']).extract(record.seq).upper()
        self.assertEqual(extracted_seq, REFERENCE.upper())

    def test_fast_suppression(self):
        record = SeqRecord(Seq(REFERENCE), id='test_ref_exact')
        entry = ENTRY_FIELDS.copy()
        
        # Create a fast hit that matches exactly
        from Bio.SeqFeature import SimpleLocation
        fast_hit_entry = entry.copy()
        fast_hit_location = SimpleLocation(0, 120, strand=1)
        fast_hits = [{'entry': fast_hit_entry, 'location': fast_hit_location}]
        
        candidates = sorted(alignment_candidates(record, [entry], fast_hits), key=lambda h: (h['coverage'],h['identity']), reverse=True)
        
        # Should not include the forward matching same locus due to suppression
        self.assertFalse(any(c['location'].strand == 1 for c in candidates))

    def test_no_zero_reverse_incidental(self):
        record = SeqRecord(Seq(REFERENCE), id='test_ref_normal')
        entry = ENTRY_FIELDS.copy()
        
        fast_hits = []
        candidates = sorted(alignment_candidates(record, [entry], fast_hits), key=lambda h: (h['coverage'],h['identity']), reverse=True)
        
        # Ensure we don't get zero reverse incidental weak hits as a failure condition
        self.assertGreater(len(candidates), 0)

