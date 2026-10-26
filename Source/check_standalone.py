"""Verify the actual distribution without installed references, Python, or Ollama."""
import hashlib
import json
import os
import shutil
import struct
import subprocess
import tempfile
from pathlib import Path
from Bio import SeqIO
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.SeqFeature import SeqFeature, SimpleLocation
from packet_writer import append_features_bytes

root = Path(__file__).resolve().parent.parent
distribution = Path(os.environ.get('PLASMID_DISTRIBUTION', str(root)))
reference = next(r for r in json.loads((distribution/'protein_features.json').read_text(encoding='utf-8'))
                 if r['canonical_name']=='URA3')
coding = ''.join(unambiguous_dna_by_id[1].back_table[a] for a in reference['sequence'])
sequence = 'TAA' + coding + 'TAA'
def packet(kind, payload):
    return bytes([kind])+struct.pack('>I',len(payload))+payload
raw = packet(9,b'SnapGene'+struct.pack('>HHH',1,1,1))+packet(0,b'\x01'+sequence.encode())
raw = append_features_bytes(raw,[SeqFeature(SimpleLocation(0,3,strand=1),type='misc_feature',qualifiers={'label':['preserved control']})])
with tempfile.TemporaryDirectory(prefix='Plasmid standalone ') as temp:
    folder = Path(temp)
    app = folder/'Portable companion α'
    app.mkdir()
    for name in ['PlasmidAssistant.exe','protein_features.json']:
        shutil.copy2(distribution/name,app/name)
    shutil.copytree(distribution/'_internal',app/'_internal')
    source, output = folder/'public reference.dna', folder/'annotated.dna'
    source.write_bytes(raw)
    env = {k:v for k,v in os.environ.items() if k not in ('PYTHONHOME','PYTHONPATH','SNAPGENE_HOME')}
    env['PATH']=''
    command=[str(app/'PlasmidAssistant.exe'),'--cli',str(source),'--output',str(output),'--accept',reference['feature_id']]
    run=subprocess.run(command,cwd=app,env=env,timeout=120,capture_output=True)
    assert run.returncode==0, f'Executable failed: {run.returncode}'
    saved=SeqIO.read(output,'snapgene')
    assert str(saved.seq)==sequence
    assert saved.features[0].qualifiers['label']==['preserved control']
    assert any(f.qualifiers.get('label')==['URA3'] for f in saved.features)
    assert source.read_bytes()==raw
    print(json.dumps({'standalone_executable':True,'installed_examples':False,'external_python':False,
        'empty_PATH':True,'protein_annotation_exported':True,'original_preserved':True,
        'exe_sha256':hashlib.sha256((app/'PlasmidAssistant.exe').read_bytes()).hexdigest()}),flush=True)
