import pathlib
import hashlib
from Bio import SeqIO
from Bio.SeqFeature import ExactPosition

CANONICAL_LABELS = [
    "AmpR", "KanR", "CmR", "HygR", "NatR", "BleR", "URA3", "HIS3", "LEU2",
    "TRP1", "ADE2", "MET15", "EGFP", "GFP", "mCherry", "mRFP1", "YFP", "CFP",
    "ori", "f1 ori", "2-micron ori", "ARS1", "CEN6", "T7 promoter", "SP6 promoter",
    "T3 promoter", "lac promoter", "tac promoter", "TEF1 promoter", "ADH1 promoter",
    "CYC1 terminator", "ADH1 terminator", "T7 terminator", "mEGFP", "mTagBFP2",
    "SYFP2", "GAP promoter", "AOX1 terminator", "PpHIS4", "p15A ori"
]

CANONICAL_LABELS += [
    'BleoR', 'EmGFP', 'Cycle 3 GFP', 'Topaz YFP', 'SV40 ori', 'ColA ori',
    'oriT', 'oriV', '2\u00b5 ori', 'CMV enhancer', 'lac operator', 'tet operator',
    'lac operator (symmetric)', 'IRES', 'RBS', 'RBS (BBa_B0034)',
    'SV40 poly(A) signal', 'HSV TK poly(A) signal', 'bGH poly(A) signal',
    'OpIE-2 poly(A) signal', 'poly(A) signal'
]

def build_library(directory, max_per_name=3):
    entries = []
    errors = []
    
    try:
        path_obj = pathlib.Path(directory)
    except Exception:
        return [], [f"Invalid directory: {directory}"]
        
    if not path_obj.exists() or not path_obj.is_dir():
        return [], [f"Path is not a directory: {directory}"]

    for path in sorted(path_obj.rglob('*.dna')):
        try:
            record = SeqIO.read(path, 'snapgene')
        except Exception as e:
            errors.append(str(path) + str(e))
            continue
            
        if not hasattr(record, 'features'):
            continue

        for feature in record.features:
            loc_str = str(feature.location).upper()
            
            # Reject missing location or CompoundLocation with fuzzy endpoints
            if feature.location is None:
                continue
                
            try:
                start = feature.location.start
                end = feature.location.end
            except Exception:
                continue

            if any(not isinstance(p.start, ExactPosition) or not isinstance(p.end, ExactPosition) for p in feature.location.parts):
                continue
                
            # Check for ambiguous sequence chars in the extracted region
            seq_str = str(feature.extract(record.seq))
            if set(seq_str.upper()) - set('ACGT'):
                continue

            if len(seq_str) < 12:
                continue
                
            # Check canonical label
            original_label = feature.qualifiers.get('label', [None])[0]
            if not original_label:
                continue
            if not (original_label.casefold() in {l.casefold() for l in CANONICAL_LABELS}
                    or 'promoter' in original_label.casefold()
                    or 'terminator' in original_label.casefold()):
                continue

            canonical_name = next((l for l in CANONICAL_LABELS if l.casefold() == original_label.casefold()), original_label)
            if not canonical_name:
                continue
                
            seq_upper = seq_str.upper()
            
            # Dedup tuple(name, sequence)
            key = (canonical_name, seq_upper)
            if any((e['canonical_name'], e['sequence']) == key for e in entries):
                continue

            # Check max_per_name limit
            current_count = sum(1 for e in entries if e['name'] == canonical_name)
            if current_count >= max_per_name:
                continue

            feature_type = getattr(feature, 'type', 'unknown')
            organism = record.annotations.get('organism', 'unknown')
            
            seq_sha256 = hashlib.sha256(seq_upper.encode()).hexdigest()[:12]
            path_sha256 = hashlib.sha256(path.read_bytes()).hexdigest()

            entry = {
                "feature_id": canonical_name + '-' + seq_sha256,
                "name": canonical_name,
                "sequence": seq_upper,
                "canonical_name": canonical_name,
                "aliases": [original_label],
                "feature_type": feature_type,
                "organism": organism,
                "expected_length": len(seq_upper),
                "description": "Extracted local reference; independent curation pending",
                "source": str(path.resolve()),
                "source_sha256": path_sha256,
                "source_location": str(feature.location),
                "version": "1"
            }
            
            entries.append(entry)

    return entries, errors
