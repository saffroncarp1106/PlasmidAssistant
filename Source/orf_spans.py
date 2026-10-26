def complete_orf_spans(sequence, circular=False, minimum=300):
    if minimum < 6: raise ValueError("minimum must be >= 6")
    seq = str(sequence).upper()
    n = len(seq)
    scan = seq * 2 if circular else seq
    stops = {"TAA", "TAG", "TGA"}
    starts = ["ATG"]
    result = []
    for frame in range(3):
        positions = range(frame, len(scan) - 2, 3)
        next_stop = None
        for pos in reversed(positions):
            codon = scan[pos:pos+3]
            if set(codon) - set('ACGT'):
                next_stop = None
                continue
            if codon in stops:
                next_stop = pos
            elif codon == starts[0] and pos < n and next_stop is not None:
                length = next_stop + 3 - pos
                if length >= minimum and (not circular or length <= n):
                    result.append((pos, length))
    return sorted(result)
