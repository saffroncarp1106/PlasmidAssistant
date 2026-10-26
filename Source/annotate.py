"""Conservative sequence import and unknown ORF review foundation."""
import argparse
import copy
import json
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from matching import load_database, match_features, to_feature, report
from export_result import save_result


from orf_detector import detect_orfs


def read_record(path):
    path = Path(path)
    fmt = {'.gb': 'genbank', '.gbk': 'genbank', '.genbank': 'genbank',
           '.fa': 'fasta', '.fasta': 'fasta', '.fna': 'fasta', '.dna': 'snapgene'}.get(path.suffix.lower())
    if fmt is None:
        raise ValueError('Supported inputs: GenBank, FASTA and SnapGene .dna.')
    if fmt == 'snapgene':
        record = SeqIO.read(path, fmt)
    else:
        with path.open(encoding='utf-8-sig') as handle:
            record = SeqIO.read(handle, fmt)
    if not record.seq or set(str(record.seq).upper()) - set('ACGTRYSWKMBDHVN'):
        raise ValueError('Expected a nonempty DNA sequence using IUPAC symbols.')
    return record


def annotate_sequence(record, minimum=300):
    result = copy.deepcopy(record)
    candidates = detect_orfs(record.seq, record.annotations.get('topology') == 'circular', minimum)
    existing_cds = [f for f in record.features if f.type == 'CDS' and f.location is not None]
    candidates = [f for f in candidates if not any(
        f.location.strand == old.location.strand and
        set(f.location).issubset(set(old.location)) for old in existing_cds)]
    return result, candidates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--circular', action='store_true', help='Set topology for FASTA input')
    parser.add_argument('--minimum-orf', type=int, default=300)
    parser.add_argument('--include-unknown-orfs', action='store_true')
    default_database = (Path(sys.executable).parent if getattr(sys, 'frozen', False) else Path(__file__).parent) / 'local_features.json'
    parser.add_argument('--database', type=Path, default=default_database, help='Reference database; bundled library loads automatically')
    parser.add_argument('--minimum-identity', type=float, default=.95)
    parser.add_argument('--fast-only', action='store_true', help='Skip partial and gapped alignment suggestions')
    parser.add_argument('--alignment-identity', type=float, default=.7)
    parser.add_argument('--alignment-coverage', type=float, default=.5)
    parser.add_argument('--max-windows', type=int, default=20)
    parser.add_argument('--accept', action='append', default=[], help='Accept reviewed feature_id')
    parser.add_argument('--reject', action='append', default=[], help='Reject feature_id')
    args = parser.parse_args()
    if args.minimum_orf < 6:
        parser.error('minimum ORF length must be at least 6 bp')
    extension = '.dna' if args.input.suffix.lower() == '.dna' else '.gb'
    output = args.output or args.input.with_name(args.input.stem + '_autoannotated' + extension)
    if output.suffix.lower() not in ('.dna', '.gb', '.gbk', '.genbank'):
        parser.error('Output must use .dna or a GenBank extension')
    input_bytes = args.input.read_bytes()
    record = read_record(args.input)
    if args.circular:
        record.annotations['topology'] = 'circular'
    if not 0 < args.minimum_identity <= 1:
        parser.error('minimum identity must be between 0 and 1')
    from analysis_service import AnnotationConfig, analyze_record
    config = AnnotationConfig(args.minimum_orf, args.minimum_identity, not args.fast_only,
                              args.alignment_identity, args.alignment_coverage, args.max_windows)
    from reference_loader import load_references
    analysis = analyze_record(record, load_references(args.database), config)
    result, candidates, hits = analysis['record'], analysis['unknown_orfs'], analysis['known_features']
    accepted = [h for h in hits if h['status'] not in ('existing', 'superseded')
                and h['entry']['feature_id'] not in args.reject
                and (h['confidence'] == 'High confidence' or h['entry']['feature_id'] in args.accept)]
    from selection_guard import validate_selection
    selected = [to_feature(h) for h in accepted]
    validate_selection(selected)
    result.features.extend(selected)
    candidates = [f for f in candidates if not any(
        f.location.strand == h['location'].strand and set(f.location).issubset(set(h['location']))
        for h in accepted)]
    if args.include_unknown_orfs:
        result.features.extend(candidates)
    result.annotations.setdefault('molecule_type', 'DNA')
    # Exclusive creation prevents overwriting inputs or previous output.
    if args.input.read_bytes() != input_bytes:
        raise ValueError('Source changed during analysis; reopen it before exporting.')
    save_result(args.input, output, record, result)
    print(json.dumps({'output': str(output), 'existing_features': len(record.features),
                      'known_features': [report(h) for h in hits],
                      'unknown_orfs': [{'location': str(f.location), 'length': len(f)}
                                       for f in candidates]}, indent=2))


if __name__ == '__main__':
    main()
