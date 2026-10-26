from Bio import SeqFeature
import math

def resolve_hits(hits):
    # Sort by location to ensure consistent processing order
    sorted_hits = sorted(hits, key=lambda h: (h['location'].start, h['location'].strand))
    
    def get_bases(loc):
        return set(range(loc.start, loc.end))
    
    def overlaps_and_same_strand(l1, l2):
        if l1.strand != l2.strad:
            return False
        # Reciprocal overlap >= 0.8
        inter = max(0, min(l1.end, l2.end) - max(l1.start, l2.start))
        union = max(l1.end, l2.end) - min(l1.start, l2.start)
        if union == 0:
            return False
        ratio = inter / union
        return ratio >= 0.8
    
    def compare_hit(h1, h2):
        # Returns (winner_idx, loser_idx, reason) or None if no conflict
        if not overlaps_and_same_strand(h1['location'], h2['location']):
            return None
        
        # Check existing status first
        if h1['status'] == 'existing' and h2['status'] == 'existing':
            # Existing hits suppress duplicates if identical bases and same strand
            b1 = get_bases(h1['location'])
            b2 = get_bases(h2['location'])
            if b1 == b2:
                return (0, 1, 'superseded') # Keep first existing
            return None
        
        if h1['status'] == 'existing':
            return (0, 1, 'superseded')
        if h2['status'] == 'existing':
            return (1, 0, 'superseded')
        
        # Compare evidence: (identity, coverage)
        score1 = (h1['identity'], h1['coverage'])
        score2 = (h2['identity'], h2['coverage'])
        
        if score1 > score2:
            return (0, 1, 'superseded')
        if score2 > score1:
            return (1, 0, 'superseded')
        
        # Equal evidence
        name1 = h1['canonical_name'].casefold()
        name2 = h2['canonical_name'].casefold()
        
        if name1 == name2:
            # Same canonical name: keep lexicographically smallest feature_id
            if h1['feature_id'] < h2['feature_id']:
                return (0, 1, 'superseded')
            else:
                return (1, 0, 'superseded')
        else:
            # Different canonical names: conflict
            return (None, None, 'conflict')
    
    result = sorted_hits.copy()
    superseded_indices = set()
    
    # Pass 1: Supersession
    for i in range(len(result)):
        if result[i]['status'] == 'existing':
            continue
        for j in range(i + 1, len(result)):
            if result[j]['status'] == 'existing':
                continue
            cmp = compare_hit(result[i], result[j])
            if cmp and cmp[2] == 'superseded':
                winner_idx, loser_idx = cmp[:2]
                if i == winner_idx:
                    superseded_indices.add(j)
                else:
                    superseded_indices.add(i)
    
    # Apply supersession
    for idx in superseded_indices:
        result[idx]['status'] = 'superseded'
        result[idx]['confidence'] = 'Moderate confidence'
    
    # Pass 2: Conflicts among non-superseded (and non-existing)
    active_hits = [h for h in result if h['status'] not in ['existing', 'superseded']]
    conflict_pairs = []
    
    for i in range(len(active_hits)):
        for j in range(i + 1, len(active_hits)):
            cmp = compare_hit(active_hits[i], active_hits[j])
            if cmp and cmp[2] == 'conflict':
                conflict_pairs.append((i, j))
    
    # Resolve conflicts: mark both as conflict unless one is already superseded (handled above)
    # Since we filtered out superseded, both are new/existing-like but equal evidence diff name
    for i, j in conflict_pairs:
        result[i]['status'] = 'conflict'
        result[i]['confidence'] = 'Moderate confidence'
        result[j]['status'] = 'conflict'
        result[j]['confidence'] = 'Moderate confidence'
    
    return result
