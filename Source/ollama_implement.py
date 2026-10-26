"""Ask the local worker for reviewable files and stage them within the project."""
import json
import urllib.request
from pathlib import Path

root = Path(__file__).resolve().parent
prompt = '''Implement a local reference database builder for this plasmid annotator.
Return JSON only: {"files": [{"path":"build_database.py","content":"..."}], "report":"..."}.
Produce build_database.py and test_database.py. Builder must read SnapGene .dna
files using Bio.SeqIO, from ../Plasmids recursively. Select ONLY explicit labels
AmpR, KanR, URA3, HIS3, LEU2, TRP1, EGFP, mCherry, ori, T7 promoter, CYC1 terminator.
Extract feature sequences in biological orientation using feature.extract.
Require unambiguous ACGT and >=12 bases; skip compound or fuzzy features for now.
Deduplicate by label and sequence; stop at 15 entries. Capture source path,
SHA256 file digest, original location, organism when available, aliases, length.
Output data/local_features.json, refusing overwrite. Required fields feature_id,
canonical_name, feature_type (original GenBank type), sequence, source (string),
version. No invented sequences. This is locally extracted reference data, not
independently curated truth. Include meaningful unittest tests using temporary
GenBank records or a testable extraction helper. Do not touch installation files.
Return assumptions and unresolved issues in report. Keep code straightforward.
'''
prompt = '''Return one JSON object with exactly these keys: files and report.
files is a list containing ONE object with path build_database.py and content
containing complete executable Python code. No markdown, no repeated keys.
Implement this exact simple algorithm using pathlib, hashlib, json and Bio.SeqIO:
main takes input directory and output JSON path as positional argparse arguments.
For sorted Path(input).rglob('*.dna'), parse SeqIO.read(path, 'snapgene').
For each record.features, label = feature.qualifiers.get('label', [''])[0].
Only labels AmpR, KanR, URA3, HIS3, LEU2, TRP1, EGFP, mCherry, ori,
T7 promoter, CYC1 terminator. Require isinstance(location, SimpleLocation)
and isinstance(location.start, ExactPosition) and likewise end. Import these
from Bio.SeqFeature. sequence = str(feature.extract(record.seq)).upper().
Skip if length <12 or set(sequence)-set('ACGT'). Dedup (label, sequence).
For each entry output feature_id as label plus first 12 sha256(sequence.encode())
hex chars, canonical_name label, feature_type feature.type, sequence,
source str(path.resolve()), source_sha256 sha256(path.read_bytes()).hexdigest(),
original_location str(feature.location), version '1', expected_length len(sequence).
Stop at 15 entries total. Output with open(output, 'x', encoding='utf-8') json.dump.
On parse error print filename and error to stderr and continue. Raise if zero entries.
Use a main guard. Do not run anything or add extra files. Keep under 90 lines.
'''
print('PRIMARY -> QWEN\n' + prompt, flush=True)
request = urllib.request.Request('http://localhost:11434/api/chat', data=json.dumps({
    'model': 'qwen3.5:9b', 'messages': [{'role': 'user', 'content': prompt}],
    'format': 'json', 'think': False, 'stream': True}).encode(),
    headers={'Content-Type': 'application/json'})
chunks = []
with urllib.request.urlopen(request, timeout=1200) as response:
    for line in response:
        content = json.loads(line).get('message', {}).get('content', '')
        chunks.append(content)
        print(content, end='', flush=True)
raw = ''.join(chunks)
(root / 'ollama-implementation.json').write_text(raw, encoding='utf-8')
result = json.loads(raw)
stage = root / 'worker-staging-v2'
stage.mkdir(exist_ok=True)
for item in result['files']:
    if item['path'] not in {'build_database.py', 'test_database.py'}:
        raise ValueError('Unexpected worker output path')
    with (stage / item['path']).open('x', encoding='utf-8') as handle:
        handle.write(item['content'])
print('\nFiles staged for primary review. ' + result.get('report', ''))
