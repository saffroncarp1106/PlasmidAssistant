def validate_selection(features):
    for i in range(len(features)):
        for j in range(i+1, len(features)):
            a, b = features[i], features[j]
            if not (a.location and b.location): continue
            if a.location.strand != b.location.strand: continue
            if a.location.strand not in (-1, 1): continue
            set_a = set(a.location)
            set_b = set(b.location)
            inter = set_a & set_b
            if not inter: continue
            if len(inter)/len(set_a) < 0.8 or len(inter)/len(set_b) < 0.8: continue
            raise ValueError('Selected annotations overlap at the same locus. Select one competing candidate before saving.')
    return None