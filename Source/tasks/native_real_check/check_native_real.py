import json, tempfile, hashlib, pathlib
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation
from dna_writer import write_annotated_dna

root = Path(__file__).resolve().parent
source = Path(json.loads((root / 'local_features.json').read_text())[0]['source'])
before = source.read_bytes()

record = SeqIO.read(source, 'snapgene')
feature = SeqFeature(SimpleLocation(10, 40, strand=-1), type='misc_feature', qualifiers={'label': ['Native writer validation'], 'note': ['Synthetic test annotation; no biological identity asserted']})

with tempfile.TemporaryDirectory() as temp:
    destination = Path(temp) / 'annotated.dna'
    write_annotated_dna(source, destination, [feature])
    output = SeqIO.read(destination, 'snapgene')
    
    assert output.seq == record.seq
    assert len(output.features) == len(record.features) + 1
    assert any(f.location == feature.location and f.qualifiers.get('label') == feature.qualifiers['label'] for f in output.features)
    
    output_bytes = destination.read_bytes()
    
    try:
        write_annotated_dna(source, destination, [feature])
    except (ValueError, FileExistsError):
        pass
    else:
        raise AssertionError("Expected exception on overwrite")
    
    assert destination.read_bytes() == output_bytes
    
    try:
        write_annotated_dna(source, source, [feature])
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError when writing to source")
    
    assert source.read_bytes() == before

print(json.dumps({
    "native_roundtrip_passed": True,
    "source_length": len(record.seq),
    "existing_features": len(record.features),
    "overwrite_refused": True
}))