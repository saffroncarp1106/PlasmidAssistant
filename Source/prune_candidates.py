def prune_candidates(hits):
    hits = sorted(hits, key=lambda h: (h['coverage'], h['identity']), reverse=True)
    kept = []
    for h in hits:
        occupied = set(h['location'])
        duplicate = any(
            h['entry']['feature_id'] == old['entry']['feature_id']
            and h['location'].strand == old['location'].strand
            and min(len(occupied), len(old['location'])) > 0
            and len(occupied & set(old['location'])) / min(len(occupied), len(old['location'])) >= .8
            for old in kept
        )
        if not duplicate:
            kept.append(h)
    return kept
