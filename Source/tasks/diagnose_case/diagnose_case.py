import copy
from matching import match_features

def diagnose(record, entries):
    names = {e['canonical_name'].casefold() for e in entries}
    expected = [f for f in record.features if f['label'].casefold() in names]
    masked = copy.deepcopy(record)
    masked.features = [f for f in record.features if f['label'].casefold() not in names]
    predicted = match_features(masked, entries)
    result = []
    for h in predicted:
        if h['status'] != 'new' or h['confidence'] != 'High confidence':
            continue
        name = h['entry']['canonical_name']
        loc = h['location']
        matches = [f for f in expected if f['label'].casefold() == name.casefold() and f['location'] == loc]
        if matches:
            continue
        overlaps = []
        for f in expected:
            if f['label'].casefold() == name.casefold():
                overlap = len(set(loc) & set(f['location'])) / max(len(loc), len(f['location']))
                overlaps.append({'location': str(f['location']), 'overlap': overlap})
        result.append({'name': name, 'location': str(loc), 'closest': overlaps})
    return result