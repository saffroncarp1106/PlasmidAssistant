"""Exercise the actual bundled executable on a local real plasmid."""
import copy
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path
from Bio import SeqIO

root = Path(__file__).resolve().parent
entries = json.loads((root / 'local_features.json').read_text())
source = Path(entries[0]['source'])
before = hashlib.sha256(source.read_bytes()).hexdigest()
record = SeqIO.read(source, 'snapgene')
exe = root / 'dist' / 'PlasmidAssistant' / 'PlasmidAssistant.exe'
with tempfile.TemporaryDirectory(prefix='plasmid portable check ') as directory:
    folder = Path(directory)
    output = folder / 'roundtrip.gb'
    subprocess.run([str(exe), '--cli', str(source), '--output', str(output)], check=True, timeout=120)
    recovered = SeqIO.read(output, 'genbank')
    assert recovered.seq == record.seq
    assert len(recovered.features) == len(record.features)
    stripped = copy.deepcopy(record)
    removed = [f for f in stripped.features if f.qualifiers.get('label') == ['ori']]
    stripped.features = [f for f in stripped.features if f not in removed]
    test_input = folder / 'missing origin.gb'
    with test_input.open('w', encoding='utf-8') as handle:
        SeqIO.write(stripped, handle, 'genbank')
    annotated = folder / 'recovered.gb'
    subprocess.run([str(exe), '--cli', str(test_input), '--output', str(annotated),
                    '--database', str(root / 'local_features.json')], check=True, timeout=180)
    result = SeqIO.read(annotated, 'genbank')
    assert removed, 'Fixture must contain an origin'
    assert all(any(f.location == target.location and f.qualifiers.get('label') == ['ori']
                   for f in result.features) for target in removed)
    assert hashlib.sha256(source.read_bytes()).hexdigest() == before
    print(json.dumps({'frozen_executable': True, 'source_length': len(record),
                      'existing_features_preserved': len(record.features),
                      'removed_origins_recovered': len(removed),
                      'original_unchanged': True,
                      'limitation': 'Reference-source smoke check, not an independent benchmark'}))
