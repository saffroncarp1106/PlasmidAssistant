import random
import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from alignment_candidates import alignment_candidates

# Fixture setup
rng = random.Random(42)
REFERENCE = ''.join(rng.choice('ACGT') for _ in range(300))

# Construct entry dict
entry = {
    'sequence': REFERENCE,
    'canonical_name': 'Synthetic variant',
    'feature_id': 'synthetic',
    'feature_type': 'CDS',
    'source': 'synthetic test',
    'version': '1'
}

# Mutate positions: range(30, 300, 4) AND range(1, 30, 4) to ensure ~25% divergence
# Preserve first 8 bases by only mutating positions >= 8
mutate_indices = set(range(30, 300, 4)) | set(range(1, 30, 4))
mutated_seq = list(REFERENCE)
for i in mutate_indices:
    if i >= 8: # Ensure seed exists
        base = REFERENCE[i]
        options = [b for b in 'ACGT' if b != base]
        mutated_seq[i] = rng.choice(options)

MUTATED_SEQ_STR = ''.join(mutated_seq)

class TestRecovery(unittest.TestCase):
    def test_75_percent_identity_recovery(self):
        record = SeqRecord(Seq(MUTATED_SEQ_STR), id='test_variant')
        
        # Call analyze_record logic via alignment_candidates directly as per source structure
        candidates = sorted(alignment_candidates(record, [entry], [], minimum_identity=.7, minimum_coverage=.8), 
                            key=lambda h: (h['coverage'], h['identity']), reverse=True)
        
        # Assert at least one hit with identity between 0.7 and 0.8
        self.assertTrue(any(0.7 <= c['identity'] < 0.8 for c in candidates))
        # Assert coverage >= 0.8
        self.assertGreaterEqual(candidates[0]['coverage'], 0.8)
        # Assert confidence is not 'High confidence'
        self.assertNotIn('High confidence', [c['confidence'] for c in candidates])

if __name__ == '__main__':
    unittest.main()