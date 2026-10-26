def prune_candidates(hits):
    hits = sorted(hits, key=lambda h: (h['coverage'], h['identity']), reverse=True)
    kept = []
    seen = set()
    for h in hits:
        loc = tuple(sorted(set(h['location']))) if isinstance(h['location'], list) else h['location']
        strand = h.get('strand', 'both')
        key = (h['feature_id'], loc, strand)
        if key not in seen:
            kept.append(h)
            seen.add(key)
    return kept