import os
import json, hashlib, subprocess, tempfile, shutil, os
from pathlib import Path
from Bio import SeqIO

root = Path(__file__).resolve().parent
entries = json.loads((root / 'local_features.json').read_text())
source = Path(entries[0]['source']).resolve()
before = source.read_bytes()
record = SeqIO.read(source, 'snapgene')
exe = Path(os.environ.get('PLASMID_ASSISTANT_EXE', str(root.parent / 'PlasmidAssistant' / 'PlasmidAssistant.exe')))

with tempfile.TemporaryDirectory(prefix='native smoke ') as folder:
    relocated = Path(folder) / 'Portable Plasmid \u03b1'
    shutil.copytree(exe.parent, relocated)
    assert not (relocated / 'local_features.json').exists()
    (Path(folder) / 'Plasmids').mkdir()
    shutil.copyfile(source, Path(folder) / 'Plasmids' / 'example.dna')
    exe = relocated / 'PlasmidAssistant.exe'
    env = os.environ.copy()
    env['PATH'] = ''
    env.pop('PYTHONHOME', None)
    env.pop('PYTHONPATH', None)
    output = Path(folder) / 'result.dna'
    subprocess.run([str(exe), '--cli', str(source), '--output', str(output), '--include-unknown-orfs'], check=True, timeout=180, env=env, cwd=relocated)
    result = SeqIO.read(output, 'snapgene')
    assert result.seq == record.seq
    assert len(result.features) > len(record.features)
    for a, b in zip(record.features, result.features):
        assert (a.type, str(a.location), a.qualifiers) == (b.type, str(b.location), b.qualifiers)
    assert source.read_bytes() == before
    print(json.dumps(dict(original_features=len(record.features), output_features=len(result.features), native_roundtrip=True, relocated=True, auto_references=True, limitation='Relocated same-machine check with empty PATH; not a clean OS test')))
