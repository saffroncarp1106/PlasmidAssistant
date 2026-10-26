from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from matching import location

def detect_orfs(sequence, circular=False, minimum=300):
    sequence = str(sequence).upper()
    n = len(sequence)
    strands = [(1, sequence), (-1, str(Seq(sequence).reverse_complement()))]
    results = []
    for strand, oriented in strands:
        orf_spans = complete_orf_spans(oriented, circular, minimum)
        for start, length in orf_spans:
            forward = start if strand == 1 else (n - start - length) % n
            loc = location(forward, length, n, strand)
            results.append(SeqFeature(loc, type='misc_feature', qualifiers={'label': ['Unknown ORF'], 'note': ['Complete ATG-to-stop ORF; identity unknown; not evidence of expression.'], 'inference': ['ab initio prediction:plasmid-assistant']}))
    return results