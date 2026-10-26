def alignment_candidates(record, entries, fast_hits, minimum_identity=0.8, minimum_coverage=0.5, max_windows=20):
    """
    Heuristic: Seeds circular scan windows, aligns best local hits, filters by identity/coverage,
    suppresses redundant fast hits (same name/strand with high overlap), and deduplicates by feature_id.
    Limitation: max_windows caps seed search depth for performance; results are partial if coverage < 0.95.
    """
    from Bio.Seq import Seq
    from matching.location import location

    sequence = record['sequence'].upper()
    n = len(sequence)
    scan = sequence * 2 if record.get('circular', False) else sequence

    results = []
    seen = set()

    for entry in entries:
        query_seq = entry['sequence']
        if len(query_seq) < 60:
            continue

        strand = entry.get('strand', 1)
        if strand == -1:
            query_seq = query_seq.reverse_complement().upper()
            strand = -1
        else:
            query_seq = query_seq.upper()

        for window in seed_windows(scan, query_seq, n, max_windows):
            left, right = window
            evidence = align_best(scan[left:right], query_seq)
            if not evidence:
                continue

            start = left + evidence['start']
            length = evidence['end'] - evidence['start']

            if not (0 < start < n and 0 < length <= n):
                continue

            identity = evidence.get('identity', 0)
            coverage = evidence.get('coverage', 0)

            if identity < minimum_identity or coverage < minimum_coverage:
                continue

            loc = location(start, length, n, strand)

            # Suppress if fast hit matches canonical name, strand, and has >=80% overlap
            for fh in fast_hits:
                if fh['entry'].get('canonical_name', '').casefold() == entry.get('canonical_name', '').casefold():
                    if fh['location'].strand == strand:
                        overlap = min(length, fh['location'].length) / max(length, fh['location'].length)
                        if overlap >= 0.8:
                            break

            key = (entry.get('feature_id'), strand, loc.str)
            if key in seen:
                continue
            seen.add(key)

            status = 'new'
            match_type = 'gapped_variant' if coverage >= 0.95 else 'partial'
            confidence = 'Moderate confidence' if identity >= 0.95 and coverage >= 0.8 else 'Weak evidence'

            results.append({
                'location': loc,
                'identity': identity,
                'coverage': coverage,
                'match': match_type,
                'confidence': confidence,
                'status': status,
                'method': 'seeded_local_alignment',
                'reference_start': start,
                'reference_end': start + length
            })

    return results