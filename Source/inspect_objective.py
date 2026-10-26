import json, Bio.SeqIO
from pathlib import Path
from analysis_service import analyze_record, AnnotationConfig
from matching import report

path = Path(r'C:\Users\llkin\Documents\WhoamI plasmid new_DAB.dna')
record = Bio.SeqIO.read(path, 'snapgene')
entries = json.loads(Path('expanded-library-optimization.json').read_text(encoding='utf-8'))
result = analyze_record(record, entries, AnnotationConfig(alignment_identity=.7))

output = {
    "existing_features": [str(f) for f in result['existing_features']],
    "known_features": [report(h) for h in result['known_features']],
    "unknown_orfs": sorted(result['unknown_orfs'], key=lambda x: len(x), reverse=True)
}

output['unknown_orfs'] = [dict(location=str(f.location), length=len(f), dna=str(f.extract(record.seq)), protein=str(f.extract(record.seq).translate(to_stop=True))) for f in output['unknown_orfs']]

with open('objective-analysis.json', 'w') as f:
    json.dump(output, f)

print(json.dumps({
    "bases": len(record.seq),
    "len_references": len(entries),
    "total_hit_count": len(result['known_features']),
    "largest_orf_lengths": [orf['length'] for orf in output["unknown_orfs"]]
}))