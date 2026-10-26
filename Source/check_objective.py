import json, hashlib, Bio.SeqIO, pathlib
from analysis_service import analyze_record
from matching import report
path = pathlib.Path(r'C:\Users\llkin\Documents\WhoamI plasmid new_DAB.dna')
before = path.read_bytes()
record = Bio.SeqIO.read(path, 'snapgene')
entries = json.loads(pathlib.Path('expanded-library-optimization.json').read_text(encoding='utf-8'))
sha256 = hashlib.sha256(before).hexdigest()
result = analyze_record(record, entries)
hits = [h for h in result['known_features'] if h['status'] != 'superseded']
summary = [report(h) for h in hits]
assert path.read_bytes() == before, 'Source bytes changed'
for name, identity, coverage in [('ASH1 (C. albicans)', .99, .99), ('mRFP1', .95, .9), ('URA3', .7, .9)]:
    assert any(h['name'] == name and h['identity'] > identity and h['coverage'] > coverage for h in summary), name
result_report = {
    "inputsha256": sha256,
    "bases": len(record.seq),
    "existing_features": [str(f) for f in result['existing_features']],
    "known_features": summary,
    "expected_feature_count": 7,
    "limitation": "Seven is assignment expectation, not forced detection; promoter identity still requires independent reference evidence"
}
pathlib.Path('objective-final-report.json').write_text(json.dumps(result_report, indent=2), encoding='utf-8')
print("Summary:", json.dumps(summary, indent=2))
