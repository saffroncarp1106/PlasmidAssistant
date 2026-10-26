from Bio.Align import PairwiseAligner

def protein_evidence(target, reference):
    if not target or not reference:
        return None
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -.5
    alignments = aligner.align(target, reference)
    alignment = next(iter(alignments), None)
    if alignment is None:
        return None
    coords = alignment.coordinates
    t0, t1 = int(coords[0, 0]), int(coords[0, -1])
    q0, q1 = int(coords[1, 0]), int(coords[1, -1])
    dt, dq = 0, 0
    denom = 0
    covered = 0
    matches = 0
    for i in range(coords.shape[1] - 1):
        a0, a1 = int(coords[0, i]), int(coords[0, i + 1])
        b0, b1 = int(coords[1, i]), int(coords[1, i + 1])
        dt, dq = a1 - a0, b1 - b0
        if dt and dq:
            if dt != dq:
                raise ValueError('Unequal diagonal')
            matches += sum(a == b and a != 'X' for a,b in zip(target[a0:a1], reference[b0:b1]))
        denom += max(dt, dq)
        covered += dq
    if denom == 0:
        identity = 0.0
    else:
        identity = matches / denom
    return {
        "start": t0,
        "end": t1,
        "reference_start": q0,
        "reference_end": q1,
        "identity": identity,
        "coverage": covered / len(reference),
        "reference_coverage": covered / len(reference),
        "query_coverage": (t1 - t0) / len(target),
        "score": float(alignment.score),
        "aligned_columns": denom
    }
