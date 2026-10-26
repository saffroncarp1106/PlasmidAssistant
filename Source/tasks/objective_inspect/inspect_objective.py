import json, Path, Bio.SeqIO
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

for orf in output["unknown_orfs"]:
    dna = str(orf.extract(record.seq))
    protein = str(orf.extract(record.seq).translate(to_stop=True))
    orf.update({"dna": dna, "protein": protein})

with open('objective-analysis.json', 'w') as f:
    json.dump(output, f)

print(json.dumps({
    "bases": len(record.seq),
    "len_references": len(entries),
    "total_hit_count": len(result['known_features']),
    "largest_orf_lengths": [len(orf.extract(record.seq)) for orf in output["unknown_orfs"]]
}))