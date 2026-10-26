"""Evidence-based feature matching with explicit reference provenance."""
import json
from resolver import resolve_hits
from seed_starts import candidate_starts
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation


def load_database(path):
    entries = json.loads(path.read_text(encoding='utf-8'))
    identifiers = set()
    for entry in entries:
        for key in ('feature_id', 'canonical_name', 'feature_type', 'sequence', 'source', 'version'):
            if not entry.get(key):
                raise ValueError(f'Missing database field: {key}')
        entry['sequence'] = entry['sequence'].upper()
        if set(entry['sequence']) - set('ACGT') or len(entry['sequence']) < 12:
            raise ValueError('References require at least 12 unambiguous DNA bases')
        if entry['feature_id'] in identifiers:
            raise ValueError('Duplicate feature_id')
        identifiers.add(entry['feature_id'])
    return entries


def location(start, length, n, strand):
    parts = [SimpleLocation(start, min(start + length, n), strand=strand)]
    if start + length > n:
        parts.append(SimpleLocation(0, start + length - n, strand=strand))
    if strand == -1:
        parts.reverse()
    return parts[0] if len(parts) == 1 else CompoundLocation(parts)


def match_features(record, entries, minimum_identity=.95):
    """Exact and ungapped full-length variants; variants require review.

    No indel or partial detection is claimed. Identity is measured over the entire
    reference. Same-locus competing names are retained for review, never guessed.
    """
    sequence = str(record.seq).upper()
    n = len(sequence)
    circular = record.annotations.get('topology') == 'circular'
    scan = sequence * 2 if circular else sequence
    hits = []
    for entry in entries:
        reference = entry['sequence']
        length = len(reference)
        if length > n:
            continue
        for strand, query in ((1, reference), (-1, str(Seq(reference).reverse_complement()))):
            for start in candidate_starts(scan, query, minimum_identity, n if circular else n - length + 1):
                matches = sum(a == b for a, b in zip(scan[start:start + length], query))
                identity = matches / length
                if identity < minimum_identity:
                    continue
                loc = location(start, length, n, strand)
                occupied = set(loc)
                existing = any(f.location is not None and f.type != 'source' and
                               f.location.strand in (strand, None, 0) and
                               occupied == set(f.location) for f in record.features)
                hits.append({'entry': entry, 'location': loc, 'identity': identity,
                             'coverage': 1.0, 'match': 'exact' if identity == 1 else 'variant',
                             'confidence': 'High confidence' if identity == 1 else 'Moderate confidence',
                             'status': 'existing' if existing else 'new'})
    return resolve_hits(hits)


def to_feature(hit):
    entry = hit['entry']
    protein = hit.get('method') in ('protein_local_alignment', 'blastx')
    units = 'amino-acid' if protein else 'DNA'
    partial = hit['match'] == 'partial' or (protein and hit['coverage'] < .99)
    return SeqFeature(hit['location'], type='misc_feature' if partial else entry['feature_type'], qualifiers={
        'label': [entry['canonical_name']],
        'note': [f"{hit['confidence']}; {hit['match']}; {units} identity={hit['identity']:.4f}; reference coverage={hit['coverage']:.4f}; not a confidence probability",
               f"Reference: {entry['feature_id']} version {entry['version']}; {entry['source']}"]
               + ([f"NCBI BLASTx; E-value={hit['evalue']}; bit score={hit['bitscore']}; standard genetic code"] if hit.get('method') == 'blastx' else [])})


def report(hit):
    return {**{k: v for k, v in hit.items() if k not in ('entry', 'location')},
            'feature_id': hit['entry']['feature_id'], 'name': hit['entry']['canonical_name'],
            'location': str(hit['location'])}
