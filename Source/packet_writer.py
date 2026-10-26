import struct, io
import xml.etree.ElementTree as ET
from Bio import SeqIO
from packet_parser import parse_packets
from feature_xml import feature_xml

def append_features_bytes(original, features):
    record = SeqIO.read(io.BytesIO(original), 'snapgene')
    packets = parse_packets(original)
    feature_packets = [p for p in packets if p[0] == 10]
    if len(feature_packets) > 1:
        raise ValueError("Multiple type-10 packets found")
    if not features:
        return original
    root = ET.fromstring(feature_packets[0][1]) if feature_packets else ET.Element('Features')
    if root.tag != 'Features':
        raise ValueError("Root tag must be 'Features'")
    next_id = 1 + max([int(e.get('id', '-1')) for e in root.findall('Feature')], default=-1)
    for f in features:
        root.append(feature_xml(f, next_id, len(record)))
        next_id += 1
    payload = ET.tostring(root, encoding='utf-8')
    new_packet = bytes([10]) + struct.pack('>I', len(payload)) + payload
    result = b''.join(new_packet if t == 10 else raw for t, payload, raw in packets)
    if not feature_packets:
        result += new_packet
    return result
