import math

def candidate_starts(scan, query, minimum_identity, max_start):
    if not query: raise ValueError("Empty query")
    if not (0 < minimum_identity <= 1): raise ValueError("Invalid identity")
    L = len(query)
    blocks = math.floor(L * (1 - minimum_identity) + 1e-9) + 1
    if blocks >= L: return [s for s in range(max_start) if s + L <= len(scan)]
    step = L // blocks
    candidates = set()
    for i in range(blocks):
        a = i * L // blocks
        b = (i + 1) * L // blocks
        block = query[a:b]
        p = scan.find(block)
        while p != -1:
            s = p - a
            if 0 <= s < max_start and s + L <= len(scan): candidates.add(s)
            p = scan.find(block, p + 1)
    return sorted(candidates)
