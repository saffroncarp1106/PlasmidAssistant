import random, json
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matching import load_database
from analysis_service import analyze_record
entries = load_database(Path(__file__).with_name('expanded-library-optimization.json'))
rng = random.Random(42017)
for i in range(10):
    sequence = ''.join((rng.choice('ACGT') for _ in range(5000)))
    record = SeqRecord(Seq(sequence), id=str(i), annotations={'topology': 'circular', 'molecule_type': 'DNA'})
    result = analyze_record(record, entries)
    automatic = [h for h in result['known_features'] if h['status'] == 'new' and h['confidence'] == 'High confidence']
    assert not automatic, 'Unexpected automatic match'
    print(json.dumps({'control':i+1,'suggestions':len(result['known_features'])}),flush=True)
print(json.dumps({'negative_controls': 10, 'bases_each': 5000, 'seed': 42017, 'automatic_calls': 0, 'limitation': 'Synthetic negatives, not biological specificity validation'}))
