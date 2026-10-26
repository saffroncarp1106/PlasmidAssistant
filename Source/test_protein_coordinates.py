import unittest, random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data.CodonTable import unambiguous_dna_by_id
from protein_candidates import protein_candidates

class CoordinateTests(unittest.TestCase):
    def setUp(self):
        self.peptide = ''.join(random.Random(5).choices(list('ACDEFGHIKLMNPQRSTVWY'),k=90))
        self.coding = ''.join(unambiguous_dna_by_id[1].back_table[a] for a in self.peptide)
        self.dna = Seq('TAA'+self.coding+'TAA')
        self.refs = [dict(feature_id='x',canonical_name='x',sequence=self.peptide)]
    def check(self,dna,strand,circular=False,wrap=False):
        record=SeqRecord(dna,annotations={'topology':'circular' if circular else 'linear'})
        hits=protein_candidates(record,[],self.refs)
        exact=[h for h in hits if h['identity']==1 and h['coverage']==1 and h['location'].strand==strand]
        self.assertTrue(exact)
        self.assertEqual(str(exact[0]['location'].extract(dna).translate(to_stop=True)),self.peptide)
        if wrap: self.assertEqual(len(exact[0]['location'].parts),2)
    def test_forward(self): self.check(self.dna,1)
    def test_reverse(self): self.check(self.dna.reverse_complement(),-1)
    def test_circular(self): self.check(self.dna[50:]+self.dna[:50],1,True,True)
    def test_circular_reverse(self): self.check((self.dna[50:]+self.dna[:50]).reverse_complement(),-1,True,True)
    def test_fusion_domain(self): self.check(Seq('GCT'*100+self.coding+'GCT'*100+'TAA'),1)

if __name__=='__main__': unittest.main()
