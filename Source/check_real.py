import json, hashlib, copy, io
from pathlib import Path
from Bio import SeqIO

root = Path(__file__).resolve().parent
with open(root / 'local_features.json') as f:
    data = json.load(f)
src_path = Path(data[0]['source'])
bytes_data = src_path.read_bytes()
sha256_before = hashlib.sha256(bytes_data).hexdigest()
record = SeqIO.read(src_path, 'snapgene')
deep_record = copy.deepcopy(record)
deep_record.annotations['molecule_type'] = 'DNA'
stream = io.StringIO()
SeqIO.write(deep_record, stream, 'genbank')
stream.seek(0)
read_back = SeqIO.read(stream, 'genbank')
assert str(read_back.seq) == str(record.seq)
assert read_back.annotations.get('topology') == record.annotations.get('topology')
assert len(read_back.features) == len(record.features)
for f1, f2 in zip(read_back.features, record.features):
    assert f1.location == f2.location
    assert f1.type == f2.type
    # GenBank represents numeric qualifier values as strings.
    normalize = lambda q: {k: [str(v) for v in values] for k, values in q.items()}
    assert normalize(f1.qualifiers) == normalize(f2.qualifiers)
assert sha256_before == hashlib.sha256(src_path.read_bytes()).hexdigest()
result = {
    'source_path': str(src_path),
    'sequence_length': len(record.seq),
    'feature_count': len(record.features),
    'roundtrip_passed': True
}
print(json.dumps(result))
