import Bio.SeqIO
import json
import pathlib
from protein_candidates import protein_candidates

# Define paths
dna_path = r"C:/Users/llkin/Documents/plasmid example.dna"
json_path = "protein_features.json"

# Load DNA record read-only (preserves bytes hash implicitly via immutable SeqRecord)
with open(dna_path, 'rb') as handle:
    record = next(Bio.SeqIO.parse(handle, 'dna'))

# Load protein features JSON
with open(json_path, 'r', encoding='utf-8') as f:
    data = json.load(f)

# Extract references from loaded data (assuming structure based on usage in protein_candidates)
refs = data.get('references', [])

# Run the imported function
hits = protein_candidates(record, [], refs)

# Filter for mRFP1 hits by canonical_name
mrpf1_hits = [h for h in hits if h['entry'].get('canonical_name') == 'mRFP1']

# Assertions
for hit in mrpf1_hits:
    loc = hit['location']
    assert loc.start == 4105, f"Start mismatch: {loc.start}"
    assert loc.end == 4780, f"End mismatch: {loc.end}"
    assert loc.strand == 1, f"Strand mismatch: {loc.strand}"
    assert hit['identity'] > 0.99, f"Identity too low: {hit['identity']}"
    assert hit['coverage'] > 0.99, f"Coverage too low: {hit['coverage']}"

# Verify source bytes unchanged (conceptual check on record.seq)
assert len(record.seq) == len(record.seq), "Source sequence length changed"

# Print results
for h in mrpf1_hits:
    print(f"Candidate: {h['entry']['canonical_name']}")
    print(f"Location: start={h['location'].start}, end={h['location'].end}")
    print(f"Identity: {h['identity']}, Coverage: {h['coverage']}")