"""Record public reference metadata for read-only sequence identification."""
import hashlib
import json
import re
from pathlib import Path
from Bio import SeqIO

wanted = {'URA3', 'HIS3', 'LEU2', 'TRP1', 'ADE2', 'MET15', 'ASH1'}
entries = []
for record in SeqIO.parse('public_references/yeast-proteome.fasta', 'fasta'):
    gene = re.search(r'\bGN=(\S+)', record.description)
    if gene and gene[1] in wanted:
        accession = record.id.split('|')[1]
        entries.append(dict(feature_id='UniProt:' + accession,
            canonical_name=gene[1], feature_type='CDS', sequence=str(record.seq),
            source='https://www.uniprot.org/uniprotkb/' + accession + '/entry',
            version='retrieved-2026-10-24', molecule_type='protein',
            description=record.description))

candidate = json.loads(Path('protein-shortlist-broad.json').read_text())[0]['candidates'][0]
assert candidate['id'] == 'sp|Q5A5Q6|ASH1_CANAL'
entries.append(dict(feature_id='UniProt:Q5A5Q6', canonical_name='ASH1 (C. albicans)',
    feature_type='CDS', sequence=candidate['sequence'],
    source='https://www.uniprot.org/uniprotkb/Q5A5Q6/entry',
    version='retrieved-2026-10-24', molecule_type='protein',
    description=candidate['description']))

record = SeqIO.read('public_references/mRFP1-CAH64892.fasta', 'fasta')
entries.append(dict(feature_id='GenBank:CAH64892.1', canonical_name='mRFP1',
    feature_type='CDS', sequence=str(record.seq),
    source='https://www.ncbi.nlm.nih.gov/protein/CAH64892.1',
    version='CAH64892.1', molecule_type='protein', description=record.description))
for entry in entries:
    entry['sequence_sha256'] = hashlib.sha256(entry['sequence'].encode()).hexdigest()
Path('protein_features.json').write_text(json.dumps(entries, indent=2), encoding='utf-8')
print([(e['canonical_name'], len(e['sequence'])) for e in entries])
