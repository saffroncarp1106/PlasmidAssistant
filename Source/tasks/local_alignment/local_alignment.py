from Bio.Align import PairwiseAligner

def align_best(target, reference):
    if not target or not reference:
        return None
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -2
    aligner.open_gap_score = -6
    aligner.extend_gap_score = -1
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
    for i in range(t1 - t0):
        a = target[t0 + i]
        b = reference[q0 + i]
        if a == b and a in "ACGT":
            matches += 1
        dt = (t0 + i + 1) - (t0 + i)
        dq = (q0 + i + 1) - (q0 + i)
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
        "score": float(alignment.score),
        "aligned_columns": denom
    }