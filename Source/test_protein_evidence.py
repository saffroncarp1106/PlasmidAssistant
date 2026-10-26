import unittest
from protein_evidence import protein_evidence

class TestProteinEvidence(unittest.TestCase):
    def test_single_insertion(self):
        reference = 'MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE'
        query = reference[:20] + 'A' + reference[20:]
        result = protein_evidence(query, reference)
        self.assertIsNotNone(result)
        self.assertAlmostEqual(result['identity'], len(reference)/(len(reference)+1), places=5)

    def test_empty_query(self):
        self.assertIsNone(protein_evidence('', 'MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE'))

    def test_empty_reference(self):
        self.assertIsNone(protein_evidence('MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE', ''))

if __name__ == '__main__':
    unittest.main()