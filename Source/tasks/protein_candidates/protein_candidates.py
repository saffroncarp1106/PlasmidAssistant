from protein_evidence import protein_evidence

def protein_candidates(record, orfs, references):
    results = []
    for orf in sorted(orfs, key=lambda x: len(x.location), reverse=True):
        protein = str(orf.extract(record.seq).translate(to_stop=True))
        for entry in references:
            evidence = protein_evidence(protein, entry['sequence'])
            if not evidence: continue
            if evidence['identity'] < 0.7 or evidence['query_coverage'] < 0.8 or evidence['reference_coverage'] < 0.7: continue
            prior = next((p for p in results if p['entry']['feature_id']==entry['feature_id'] and p['location'].strand==orf.location.strand and p['location'].issubset(orf.location)), None)
            if prior: continue
            results.append({**entry, 'location': orf.location, **evidence, 'match': 'protein_variant', 'confidence': 'Moderate confidence' if evidence['identity']>=0.9 else 'Weak evidence', 'status': 'new', 'method': 'protein_local_alignment'})
    return results