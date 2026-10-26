"""One live integration check using a public UniProt reference, not user DNA."""
import json
from pathlib import Path
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from blast_transport import BlastClient
from blast_parser import blast_candidates

root = Path(__file__).resolve().parent
entry = next(e for e in json.loads((root/'protein_features.json').read_text(encoding='utf-8'))
             if e['canonical_name']=='URA3')
table = CodonTable.unambiguous_dna_by_id[1].back_table
sequence = ''.join(table[aa] for aa in entry['sequence'])
print('Submitting back-translation of public UniProt P03962 for BLASTx integration check.', flush=True)
rid, xml = BlastClient().search(sequence)
(root/'blast-live-result.xml').write_text(xml, encoding='utf-8')
hits = blast_candidates(xml, SeqRecord(Seq(sequence)), rid)
assert any(h['identity'] >= .99 and h['coverage'] >= .99 for h in hits), 'Expected reference match absent'
report = {'rid':rid,'public_reference':entry['source'],'user_sequence_uploaded':False,
          'candidates':len(hits),'best_identity':max(h['identity'] for h in hits)}
(root/'blast-live-report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report),flush=True)
