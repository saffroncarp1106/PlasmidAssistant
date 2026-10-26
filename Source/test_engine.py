import io
import unittest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation
from annotate import detect_orfs, annotate_sequence
from matching import match_features, to_feature

REFERENCE = 'ATGACCGTACGATTCAGGTAA'
def entry(name='synthetic'):
    return dict(feature_id=name, canonical_name=name, feature_type='misc_feature', sequence=REFERENCE, source='Synthetic test fixture', version='1')
def record(seq, circular=False):
    return SeqRecord(Seq(seq), id='fixture', annotations={'molecule_type':'DNA', 'topology':'circular' if circular else 'linear'})

class OrfTests(unittest.TestCase):
    def test_forward_reverse(self):
        for strand, seq in [(1, 'CCCATGAAATAACCC'), (-1, str(Seq('CCCATGAAATAACCC').reverse_complement()))]:
            hits = detect_orfs(seq, minimum=9)
            self.assertEqual(len(hits), 1)
            self.assertEqual(hits[0].location.strand, strand)
            self.assertEqual(str(hits[0].extract(Seq(seq))), 'ATGAAATAA')
    def test_circular_both_strands(self):
        for seq in ['AAATAACCCATG', str(Seq('AAATAACCCATG').reverse_complement())]:
            hits = detect_orfs(seq, circular=True, minimum=9)
            self.assertEqual(len(hits), 1)
            self.assertEqual(str(hits[0].extract(Seq(seq))), 'ATGAAATAA')
            self.assertEqual(len(hits[0].location.parts), 2)
    def test_ambiguity_breaks(self):
        self.assertEqual(detect_orfs('ATGNNNAAATAA', minimum=6), [])
    def test_nested(self):
        hits = detect_orfs('ATGATGAAATAA', minimum=6)
        self.assertEqual(sorted(len(h) for h in hits), [9, 12])
    def test_existing_and_deep_copy(self):
        rec = record('ATGAAATAA')
        rec.features = [SeqFeature(SimpleLocation(0,9,strand=1), type='CDS', qualifiers={'label':['original']})]
        result, candidates = annotate_sequence(rec, minimum=9)
        self.assertEqual(candidates, [])
        result.features[0].qualifiers['label'][0] = 'changed'
        self.assertEqual(rec.features[0].qualifiers['label'], ['original'])
    def test_genbank_roundtrip(self):
        rec = record('AAATAACCCATG', True)
        rec.features = detect_orfs(rec.seq, True, 9)
        stream = io.StringIO()
        SeqIO.write(rec, stream, 'genbank')
        stream.seek(0)
        recovered = SeqIO.read(stream, 'genbank')
        self.assertEqual(recovered.seq, rec.seq)
        self.assertEqual(recovered.features[0].location, rec.features[0].location)
        self.assertEqual(recovered.features[0].qualifiers, rec.features[0].qualifiers)

class MatchingTests(unittest.TestCase):
    def test_exact_and_reverse(self):
        for strand, seq in [(1, REFERENCE), (-1, str(Seq(REFERENCE).reverse_complement()))]:
            hits = match_features(record(seq), [entry()], 1)
            self.assertEqual(len(hits), 1)
            self.assertEqual(hits[0]['location'].strand, strand)
            self.assertEqual(hits[0]['confidence'], 'High confidence')
            self.assertEqual(str(to_feature(hits[0]).extract(Seq(seq))), REFERENCE)
    def test_circular(self):
        seq = REFERENCE[7:] + 'CCCC' + REFERENCE[:7]
        hits = match_features(record(seq, True), [entry()], 1)
        self.assertEqual(len(hits), 1)
        self.assertEqual(str(to_feature(hits[0]).extract(Seq(seq))), REFERENCE)
        self.assertEqual(match_features(record(seq), [entry()], 1), [])
    def test_variant_requires_review(self):
        seq = REFERENCE[:5] + 'A' + REFERENCE[6:]
        hit, = match_features(record(seq), [entry()], .95)
        self.assertEqual(hit['match'], 'variant')
        self.assertEqual(hit['confidence'], 'Moderate confidence')
    def test_competing_references(self):
        hits = match_features(record(REFERENCE), [entry('a'), entry('b')], 1)
        self.assertEqual(len(hits), 2)
        self.assertTrue(all(h['status']=='conflict' for h in hits))
    def test_existing(self):
        rec = record(REFERENCE)
        rec.features = [SeqFeature(SimpleLocation(0,len(REFERENCE),strand=1), type='misc_feature')]
        hit, = match_features(rec, [entry()], 1)
        self.assertEqual(hit['status'], 'existing')

if __name__ == '__main__':
    unittest.main()
