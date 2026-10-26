import unittest
from Bio.Seq import Seq
from translated_segments import translated_segments

class SegmentTests(unittest.TestCase):
    def test_missing_atg_and_reverse(self):
        sequence = 'CCCGCTGCTGCTTAA'
        for dna, strand in [(sequence, 1), (str(Seq(sequence).reverse_complement()), -1)]:
            hits = translated_segments(dna, minimum_aa=3)
            self.assertTrue(any(h['start']==0 and h['strand']==strand and h['protein']=='PAAA' and h['has_stop'] for h in hits))
    def test_boundaries(self):
        for boundary in ['TAA', 'NNN']:
            hits = translated_segments('GCTGCT'+boundary+'GCTGCTTAA', minimum_aa=2)
            for start in [0,9]:
                self.assertTrue(any(h['start']==start and h['strand']==1 and h['protein']=='AA' for h in hits))
            self.assertTrue(all('*' not in h['protein'] and 'X' not in h['protein'] for h in hits))
    def test_circular_size(self):
        dna='GCTGCTTAAGCTGCTTAA'
        self.assertTrue(all(len(h['protein'])*3<=len(dna) for h in translated_segments(dna,True,2)))

if __name__=='__main__': unittest.main()
