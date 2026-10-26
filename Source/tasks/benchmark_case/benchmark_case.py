from copy import deepcopy
from collections import Counter
import matching

def evaluate_record(record, entries):
    names = {e['canonical_name'].casefold() for e in entries}
    indices = [i for i, f in enumerate(record.features) if f.location and f.qualifiers.get('label', [''])[0].casefold() in names]
    if not indices: return None
    masked = deepcopy(record)
    masked.features = [f for i, f in enumerate(masked.features) if i not in indices]
    hits = matching.match_features(masked, entries)
    def key(name, loc): return (name.casefold(), loc.strand, tuple(sorted(set(loc))))
    expected = Counter(key(record.features[i].qualifiers['label'][0], record.features[i].location) for i in indices)
    actual = Counter(key(h['entry']['canonical_name'], h['location']) for h in hits if h['status'] == 'new' and h['confidence'] == 'High confidence')
    return {'removed': sum(expected.values()), 'tp': sum((expected & actual).values()), 'fp': sum((actual - expected).values()), 'fn': sum((expected - actual).values())}