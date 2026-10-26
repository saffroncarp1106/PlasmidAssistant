import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from analysis_service import analyze_record, AnnotationConfig
from test_alignment_candidates import REFERENCE, ENTRY_FIELDS

class TestAnalyzeRecord(unittest.TestCase):
    def test_partial_input(self):
        record = SeqRecord(Seq(REFERENCE[20:110]), id='partial')
        result = analyze_record(record, [ENTRY_FIELDS])
        self.assertTrue(any(h['match'] == 'partial' for h in result['known_features']))
        self.assertGreaterEqual(result['known_features'][0]['coverage'], 0.7)
        self.assertNotIn('High confidence', [h['confidence'] for h in result['known_features']])
        self.assertEqual(len(result['existing_features']), 0)

    def test_negative_record(self):
        record = SeqRecord(Seq('N'*1000), id='noise')
        result = analyze_record(record, [ENTRY_FIELDS])
        self.assertEqual(result['known_features'], [])
        self.assertEqual(result['unknown_orfs'], [])

    def test_invalid_config(self):
        with self.assertRaises(ValueError):
            config = AnnotationConfig(alignment_identity=1.1)
            analyze_record(SeqRecord(Seq('A')), [], config=config)

if __name__ == '__main__':
    unittest.main()