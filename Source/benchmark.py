"""Synthetic correctness benchmark; does not estimate real-world accuracy."""
import json
from time import perf_counter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matching import match_features


def main():
    ref = 'GATCTACCGTATGCCAGTACCGATGCTAGC'
    entries = [dict(feature_id='synthetic-1', canonical_name='Synthetic fixture', feature_type='misc_feature', sequence=ref, source='Artificial benchmark sequence; no biological identity asserted', version='1')]
    cases = [
        ('forward', 'CCCC'+ref+'AAAA', False, {(4,1)}),
        ('reverse', 'CCCC'+str(Seq(ref).reverse_complement())+'AAAA', False, {(4,-1)}),
        ('circular', ref[10:]+'AAAA'+ref[:10], True, {(len(ref)-10+4,1)}),
        ('negative', 'A'*100, False, set()),
        ('ambiguous', ref[:10]+'N'+ref[11:], False, set()),
    ]
    results=[]
    started=perf_counter()
    for name, sequence, circular, expected in cases:
        rec=SeqRecord(Seq(sequence), annotations={'topology':'circular' if circular else 'linear'})
        hits=match_features(rec, entries, 1)
        actual={(int(h['location'].parts[0].start), h['location'].strand) for h in hits}
        results.append(dict(case=name, passed=actual==expected, true_positive=len(actual & expected), false_positive=len(actual-expected), false_negative=len(expected-actual)))
    print(json.dumps(dict(validation='SYNTHETIC ONLY: not real-world biological validation', elapsed_seconds=perf_counter()-started, cases=results), indent=2))
    return 0 if all(r['passed'] for r in results) else 1

if __name__=='__main__':
    raise SystemExit(main())
