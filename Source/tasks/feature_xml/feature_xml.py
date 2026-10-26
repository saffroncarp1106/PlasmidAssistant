from Bio.SeqFeature import SeqFeature, ExactPosition
import xml.etree.ElementTree as ET

def feature_xml(feature, feature_id, sequence_length):
    if feature.location is None:
        raise ValueError("Location cannot be None")
    if feature.location.strand not in (1, -1):
        raise ValueError(f"Strand must be 1 or -1, got {feature.location.strand}")
    
    # Check all parts have same strand and are ExactPosition
    parts = list(feature.location.parts)
    if len(parts) == 0:
        raise ValueError("No location parts found")
    
    first_strand = parts[0].strand
    for part in parts:
        if part.strand != first_strand:
            raise ValueError(f"All parts must have the same strand, got {part.strand}")
        if not isinstance(part.start, ExactPosition) or not isinstance(part.end, ExactPosition):
            raise ValueError("Both endpoints must be ExactPosition")
        if part.start >= sequence_length or part.end > sequence_length:
            raise ValueError(f"Endpoints exceed sequence_length {sequence_length}")
        if part.start < 0 or part.end <= part.start:
            raise ValueError(f"Invalid range: start={part.start}, end={part.end}")
    
    # Determine directionality
    directionality = '2' if feature.location.strand == -1 else '1'
    
    # Build Feature element
    feat_elem = ET.Element('Feature')
    feat_elem.set('id', str(feature_id))
    label = feature.qualifiers.get('label', [feature.type])[0]
    feat_elem.set('name', label)
    feat_elem.set('type', feature.type)
    feat_elem.set('directionality', directionality)
    
    # Process qualifiers
    for key, value in feature.qualifiers.items():
        if key == 'label':
            continue  # Already handled as 'name'
        if isinstance(value, list):
            val_list = value
        else:
            val_list = [value]
        for v in val_list:
            q_elem = ET.SubElement(feat_elem, 'Q')
            q_elem.set('name', str(key))
            q_elem.text = str(v)
    
    # Build Segment nodes
    segments = []
    for part in parts:
        seg_range = f'{int(part.start)+1}-{int(part.end)}'
        seg_elem = ET.SubElement(feat_elem, 'Segment')
        seg_elem.set('range', seg_range)
        seg_elem.set('type', 'standard')
        seg_elem.set('color', '#4CAF50')
    
    # Reverse segments if strand is -1
    if feature.location.strand == -1:
        segments.reverse()
    
    return feat_elem