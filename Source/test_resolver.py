import unittest
from Bio.SeqFeature import SimpleLocation, CompoundLocation
from resolver import resolve_hits

def hit(name, identity=1, location=None, status='new', canonical=None):
    return dict(entry=dict(feature_id=name, canonical_name=canonical or name),
                identity=identity, coverage=1, status=status, confidence='High confidence',
                location=location if location is not None else SimpleLocation(0, 100, strand=1))

class ResolverTests(unittest.TestCase):
    def test_same_name_different_boundaries_require_review(self):
        for reverse_order in (False, True):
            a = hit('a', canonical='ori')
            b = hit('b', canonical='ori', location=SimpleLocation(0, 95, strand=1))
            resolve_hits([b, a] if reverse_order else [a, b])
            self.assertEqual([a['status'], b['status']], ['conflict', 'conflict'])
            self.assertEqual(a['confidence'], 'Moderate confidence')

    def test_exact_wins(self):
        a, b = hit('a', .98), hit('b')
        resolve_hits([a,b])
        self.assertEqual(a['status'], 'superseded')
        self.assertEqual(b['status'], 'new')

    def test_conflict_after_suppression(self):
        hits = [hit('weak', .96), hit('one'), hit('two')]
        resolve_hits(hits)
        self.assertEqual([h['status'] for h in hits], ['superseded', 'conflict', 'conflict'])

    def test_circular_and_existing(self):
        loc = CompoundLocation([SimpleLocation(90,100,strand=1), SimpleLocation(0,10,strand=1)])
        a, b = hit('a', location=loc, status='existing'), hit('b', location=loc)
        resolve_hits([a,b])
        self.assertEqual(a['status'], 'existing')
        self.assertEqual(b['status'], 'superseded')

    def test_same_name_and_opposite_strand(self):
        a, b = hit('b',canonical='AmpR'),hit('a',canonical='ampr')
        c = hit('other',location=SimpleLocation(0,100,strand=-1))
        resolve_hits([a,b,c])
        self.assertEqual([h['status'] for h in [a,b,c]], ['superseded','new','new'])
