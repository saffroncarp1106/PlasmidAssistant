def seed_windows(scan, query, molecule_length, max_windows=20):
    windows = []
    offset = 0
    margin = max(30, len(query) // 5)
    while offset + 12 <= len(query):
        seed = query[offset:offset+12]
        pos = scan.find(seed, 0)
        while pos != -1:
            start = pos - offset
            left = max(0, start - margin)
            right = min(len(scan), start + len(query) + margin)
            if left < right and left < molecule_length:
                window = (left, right)
                if window not in windows:
                    windows.append(window)
                    if len(windows) >= max_windows:
                        return sorted(windows)
            pos = scan.find(seed, pos + 1)
        offset += 12
    return sorted(windows)