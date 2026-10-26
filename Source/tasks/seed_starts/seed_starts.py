import math

def candidate_starts(scan, query, minimum_identity, max_start):
    if not query: raise ValueError("Empty query")
    if not (0 < minimum_identity <= 1): raise ValueError("Invalid identity")
    L = len(query)
    blocks = int(math.floor(L * (1 - minimum_identity)) + 1e-9) + 1
    if blocks >= L: return [s for s in range(max_start) if s + L <= len(scan)]
    step = L // blocks
    candidates = set()
    for i in range(blocks):
        a = i * step
        b = min((i + 1) * step, L)
        block = query[a:b]
        p = scan.find(block, max_start if i == blocks - 1 else 0)
        while p != -1:
            s = p - a
            if 0 <= s and s + L <= len(scan): candidates.add(s)
            p = scan.find(block, p + 1)
    return sorted(candidates)