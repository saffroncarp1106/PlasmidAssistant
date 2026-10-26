import argparse, json, sys
from pathlib import Path
from hashlib import sha256
from Bio import SeqIO
from Bio.SeqFeature import SimpleLocation, ExactPosition

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_dir', type=Path)
    parser.add_argument('output_path', type=Path)
    args = parser.parse_args()
    
    valid_labels = {'AmpR', 'KanR', 'URA3', 'HIS3', 'LEU2', 'TRP1', 'EGFP', 'mCherry', 'ori', 'T7 promoter', 'CYC1 terminator'}
    entries = []
    
    for path in sorted(args.input_dir.rglob('*.dna')):
        if len(entries) >= 15:
            break
        try:
            record = SeqIO.read(path, 'snapgene')
        except Exception as e:
            print(f"Error parsing {path}: {e}", file=sys.stderr)
            continue
        
        for feature in record.features:
            if not isinstance(feature.location, SimpleLocation):
                continue
            loc = feature.location
            if not (isinstance(loc.start, ExactPosition) and isinstance(loc.end, ExactPosition)):
                continue
            try:
                label = feature.qualifiers.get('label', [''])[0]
                if label not in valid_labels:
                    continue
                sequence = str(feature.extract(record.seq)).upper()
                if len(sequence) < 12 or set(sequence) - {'A','C','G','T'}:
                    continue
            except Exception:
                continue
            
            if len(entries) >= 15:
                break
            key = (label, sequence)
            if key in [(e['canonical_name'], e['sequence']) for e in entries]:
                continue
            
            feature_id = label + sha256(sequence.encode()).hexdigest()[:12]
            canonical_name = label
            source_str = str(path.resolve())
            source_sha = sha256(path.read_bytes()).hexdigest()
            
            entries.append({
                'feature_id': feature_id,
                'canonical_name': canonical_name,
                'feature_type': feature.type,
                'sequence': sequence,
                'source': source_str,
                'source_sha256': source_sha,
                'original_location': str(feature.location),
                'version': '1',
                'expected_length': len(sequence)
            })
    
    if not entries:
        raise ValueError("No valid entries found.")
    
    with open(args.output_path, 'x', encoding='utf-8') as f:
        json.dump(entries, f)

if __name__ == '__main__':
    main()
