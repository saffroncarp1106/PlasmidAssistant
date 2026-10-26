from Bio.Seq import Seq
from local_alignment import align_best
import matching.location

def alignment_candidates(record, entries, fast_hits, minimum_identity=.8, minimum_coverage=.5, max_windows=20):
    """
    Finds alignment candidates using a seed heuristic.
    Note: Seed heuristic may miss divergent features due to fixed window seeds.
    """
    scan = str(record.seq).upper()
    n = len(scan)
    circular = record.annotations.get('topology') == 'circular'
    if circular:
        scan *= 2
    
    hits = []
    
    for entry in entries:
        ref_seq = entry['seq']
        if len(ref_seq) < 60:
            continue
            
        query = ref_seq
        strand_options = [(1, query), (-1, str(Seq(query).reverse_complement()).upper())]
        
        windows = set()
        for strand, q_seq in strand_options:
            q_len = len(q_seq)
            # Seed search with 12-mer offset step
            for offset in range(0, n - 11, 12):
                seed = scan[offset:offset+12]
                if not seed: continue
                for pos in range(offset, n - q_len + 1):
                    if scan[pos:pos+12] == seed:
                        predicted = pos - offset
                        margin = max(30, q_len // 5)
                        left = max(0, predicted - margin)
                        right = min(n, predicted + q_len + margin)
                        if left < right and left < n:
                            windows.add((left, right))
                        
                        if len(windows) >= max_windows:
                            break
                else:
                    continue
                break
        
        for left, right in sorted(windows):
            evidence = align_best(scan[left:right], query)
            if not evidence:
                continue
            
            start = left + evidence['start']
            length = evidence['end'] - evidence['start']
            
            if not (start < n and 0 < length <= n):
                continue
                
            identity = evidence.get('identity', 0.0)
            coverage = evidence.get('coverage', 0.0)
            
            if identity < minimum_identity or coverage < minimum_coverage:
                continue
            
            loc = matching.location.Location(start, length, n, strand)
            
            # Check fast_hits for overlap/conflict
            skip = False
            for fh in fast_hits:
                if fh.get('canonical_name', '').casefold() == entry.get('canonical_name', '').casefold():
                    if fh['strand'] == strand and fh['overlap_ratio'] >= 0.8:
                        skip = True
                        break
            
            if skip:
                continue
                
            # Dedup key
            dedup_key = (entry.get('feature_id'), strand, tuple((int(p.start), int(p.end)) for p in loc.parts))
            
            hit_dict = {
                'entry': entry,
                'location': loc,
                'identity': identity,
                'coverage': coverage,
                'match': 'partial' if coverage < 0.95 else 'gapped_variant',
                'confidence': 'Moderate confidence' if identity >= 0.95 and coverage >= 0.8 else 'Weak evidence',
                'status': 'new',
                'method': 'seeded_local_alignment',
                'reference_start': start,
                'reference_end': start + length
            }
            
            hits.append(hit_dict)
            
    return hits