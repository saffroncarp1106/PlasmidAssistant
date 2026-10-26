"""Conservative resolution of overlapping reference matches."""
def resolve_hits(hits):
    bases = [set(h['location']) for h in hits]
    scores = [(h['identity'], h['coverage']) for h in hits]

    def competing(i, j):
        if i == j or hits[i]['location'].strand != hits[j]['location'].strand:
            return False
        overlap = len(bases[i] & bases[j])
        return bool(bases[i] and bases[j]) and overlap / len(bases[i]) >= .8 and overlap / len(bases[j]) >= .8

    for i, hit in enumerate(hits):
        if hit['status'] == 'existing':
            continue
        for j, other in enumerate(hits):
            if not competing(i, j):
                continue
            same_name = hit['entry']['canonical_name'].casefold() == other['entry']['canonical_name'].casefold()
            duplicate_existing = other['status'] == 'existing' and bases[i] == bases[j]
            better = scores[j] > scores[i]
            # IDs can deduplicate identical loci, not choose biological boundaries.
            tie = scores[j] == scores[i] and same_name and bases[i] == bases[j] and other['entry']['feature_id'] < hit['entry']['feature_id']
            if duplicate_existing or better or tie:
                hit['status'] = 'superseded'
                hit['confidence'] = 'Moderate confidence'
                break
    for i, hit in enumerate(hits):
        if hit['status'] in ('existing', 'superseded'):
            continue
        for j, other in enumerate(hits):
            if other['status'] in ('existing', 'superseded') or not competing(i, j):
                continue
            if scores[i] == scores[j] and (bases[i] != bases[j] or hit['entry']['canonical_name'].casefold() != other['entry']['canonical_name'].casefold()):
                hit['status'] = 'conflict'
                hit['confidence'] = 'Moderate confidence'
    return hits
