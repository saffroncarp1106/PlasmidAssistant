from Bio.Align import PairwiseAligner

def protein_evidence(query, reference):
    if not query or not reference: return None
    aligner = PairwiseAligner(mode="local", match_score=2, mismatch_score=-1, open_gap_score=-5, extend_gap_score=-0.5)
    alignment = next(iter(aligner.align(query, reference)), None)
    if not alignment: return None
    coords = list(alignment.coordinates)
    if not coords: return None
    last_q, last_r = coords[-1]
    first_q, first_r = coords[0]
    delta_q = len(query) - (last_q - first_q + 1)
    delta_r = len(reference) - (last_r - first_r + 1)
    matches = sum(1 for k in range(last_q - first_q + 1) if query[first_q + k] == reference[first_r + k] and query[first_q + k] != 'X')
    denominator = max(delta_q, delta_r)
    return {
        "identity": matches / denominator if denominator else 0.0,
        "query_coverage": (last_q - first_q + 1) / len(query),
        "reference_coverage": (last_r - first_r + 1) / len(reference),
        "query_start": first_q, "query_end": last_q,
        "reference_start": first_r, "reference_end": last_r,
        "score": float(alignment.score)
    }