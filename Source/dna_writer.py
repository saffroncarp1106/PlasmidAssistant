import pathlib
import io
import json
from collections import Counter
from Bio import SeqIO
from packet_writer import append_features_bytes, parse_packets

def write_annotated_dna(source, destination, features):
    source = pathlib.Path(source)
    destination = pathlib.Path(destination)
    
    if destination.suffix.lower() != '.dna':
        raise ValueError("Destination must have .dna suffix")
    if source.resolve() == destination.resolve():
        raise ValueError("Source and destination cannot be identical")
    if destination.exists():
        raise ValueError("Destination already exists")
    
    original_bytes = source.read_bytes()
    before = SeqIO.read(io.BytesIO(original_bytes), 'snapgene')
    proposed = append_features_bytes(original_bytes, features)
    after = SeqIO.read(io.BytesIO(proposed), 'snapgene')
    
    if str(before.seq) != str(after.seq):
        raise ValueError("Sequence mismatch")
    if before.annotations.get('topology') != after.annotations.get('topology'):
        raise ValueError("Topology mismatch")
    
    def get_sig(f):
        return (f.type, str(f.location), json.dumps(f.qualifiers, sort_keys=True))
    
    original_sigs = Counter(map(get_sig, before.features))
    after_sigs = Counter(map(get_sig, after.features))
    
    if not original_sigs <= after_sigs:
        raise ValueError("Original features not subset of result")
    
    if len(after.features) != len(before.features) + len(features):
        raise ValueError("Feature count mismatch")
    
    feature_triples = [(f.type, str(f.location), tuple(f.qualifiers.get('label', [f.type]))) for f in features]
    proposed_triples = Counter(map(lambda x: (x[0], str(x[1]), tuple(x[2])), feature_triples))
    
    def triple(f):
        return (f.type, str(f.location), tuple(f.qualifiers.get('label', [f.type])))
    added_triples = Counter(map(triple, after.features)) - Counter(map(triple, before.features))
    if not proposed_triples <= added_triples:
        raise ValueError("Supplied features not found in result")
    
    raw_packets_orig = [p[2] for p in parse_packets(original_bytes) if p[0] != 10]
    raw_packets_prop = [p[2] for p in parse_packets(proposed) if p[0] != 10]
    
    if raw_packets_orig != raw_packets_prop:
        raise ValueError("Raw packet comparison failed")
    
    original_bytes_check = source.read_bytes()
    if original_bytes_check != original_bytes:
        raise ValueError("Source bytes changed during read")
    
    with destination.open('xb') as handle:
        handle.write(proposed)
