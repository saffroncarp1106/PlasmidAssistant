def seed_windows(scan, query, molecule_length, max_windows=20, seed_size=12, step=None):
    if not (max_windows >= 1 and seed_size >= 1 and (step is None or step >= 1)): raise ValueError("Invalid parameters")
    if step is None: step = seed_size
    margin = max(30, len(query) // 5)
    groups = {}
    for offset in range(0, len(query) - seed_size + 1, step):
        seed = query[offset:offset+seed_size]
        if scan.count(seed) > 64: continue
        pos = scan.find(seed, 0)
        while pos != -1:
            start = pos - offset
            key = start // margin
            groups.setdefault(key, []).append(start)
            pos = scan.find(seed, pos + 1)
    ranked_keys = sorted(groups.keys(), key=lambda k: (-len(groups[k]), groups[k][0]))
    windows = []
    for key in ranked_keys:
        if len(windows) >= max_windows: break
        center = groups[key][0]
        left = max(0, center - margin)
        right = min(len(scan), center + len(query) + margin)
        if left < right and left < molecule_length:
            window = (left, right)
            if window not in windows: windows.append(window)
    return sorted(windows)