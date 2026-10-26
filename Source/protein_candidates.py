from translated_segments import translated_segments
from protein_evidence import protein_evidence
from matching import location
from prune_candidates import prune_candidates

def protein_candidates(record, orfs, references):
    results = []
    
    segments = translated_segments(
        record.seq, 
        circular=record.annotations.get("topology")=="circular",
        minimum_aa=60
    )
    
    for segment in segments:
        protein = segment["protein"]
        
        for entry in references:
            evidence = protein_evidence(protein, entry['sequence'])
            
            if not evidence:
                continue
            
            if evidence['identity'] < 0.7 or evidence['end']-evidence['start'] < 60:
                continue
            if evidence['reference_coverage'] < 0.7:
                continue
            
            oriented_start = segment["start"] + 3 * evidence["start"]
            
            length = 3 * (evidence["end"] - evidence["start"])
            
            if evidence["end"] == len(protein) and segment["has_stop"]:
                length += 3
            
            n = len(record.seq)
            
            if not (oriented_start < n and 0 < length <= n):
                continue
            
            strand = segment.get("strand", 1) # Default to forward if unknown
            
            if strand == 1:
                forward_start = oriented_start
            else:
                forward_start = (n - oriented_start - length) % n
            
            location_strand = 1 if strand == 1 else -1
            
            loc_obj = location(forward_start, length, n, strand)

            result_entry = {
                'entry': entry, 
                'location': loc_obj, 
                **evidence, 
                'match': 'protein_variant', 
                'confidence': 'Moderate confidence' if evidence['identity']>=0.9 else 'Weak evidence', 
                'status': 'new', 
                'method': 'protein_local_alignment'
            }
            
            results.append(result_entry)

    return prune_candidates(results)
