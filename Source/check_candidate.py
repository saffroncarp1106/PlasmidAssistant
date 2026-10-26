import json, hashlib, subprocess, tempfile
from pathlib import Path
from Bio import SeqIO

root = Path(__file__).resolve().parent
entries = json.loads((root / 'local_features.json').read_text())
source = Path(entries[0]['source'])
before = source.read_bytes()
record = SeqIO.read(source, 'snapgene')
exe = root / 'release-candidate' / 'PlasmidAssistant' / 'PlasmidAssistant.exe'

with tempfile.TemporaryDirectory(prefix='native smoke ') as folder:
    output = Path(folder) / 'result.dna'
    subprocess.run([str(exe), '--cli', str(source), '--output', str(output), '--include-unknown-orfs'], check=True, timeout=180)
    result = SeqIO.read(output, 'snapgene')
    assert result.seq == record.seq
    assert len(result.features) >= len(record.features)
    for a, b in zip(record.features, result.features):
        assert (a.type, str(a.location), a.qualifiers) == (b.type, str(b.location), b.qualifiers)
    assert source.read_bytes() == before
    print(json.dumps(dict(original_features=len(record.features), output_features=len(result.features), native_roundtrip=True, limitation='Same-machine smoke check only')))