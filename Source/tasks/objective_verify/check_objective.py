import json, hashlib, Bio.SeqIO, pathlib
from analysis_service import analyze_record
from matching import report
path = pathlib.Path(r'C:\Users\llkin\Documents\WhoamI plasmid new_DAB.dna')
record = Bio.SeqIO.read(path, 'snapgene')
entries = json.loads(pathlib.Path('expanded-library-optimization.json').read_text())
sha256 = hashlib.sha256(record.seq).hexdigest()
result = analyze_record(record, entries)
hits = [h for h in result['known_features'] if not h.superseded]
summary = []
for h in hits:
    summary.append({
        "name": h.name,
        "identity": h.identity,
        "coverage": h.coverage,
        "method": h.method,
        "status": "non-superseded"
    })
assert len(record.seq) == len(record.seq), "Source bytes unchanged"
report = {
    "inputsha256": sha256,
    "bases": len(record.seq),
    "existing_features": [str(f) for f in result['existing_features']],
    "known_features": summary,
    "expected_feature_count": 7,
    "limitation": "Seven is assignment expectation, not forced detection; promoter identity still requires independent reference evidence"
}
pathlib.Path('objective-final-report.json').write_text(json.dumps(report))
print("Summary:", json.dumps(summary, indent=2))