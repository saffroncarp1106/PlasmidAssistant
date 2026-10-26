import json
from pathlib import Path
from Bio import SeqIO
from protein_evidence import protein_evidence

def main():
    obj = json.loads(Path("objective-analysis.json").read_text())
    shortlist = json.loads(Path("protein-shortlist.json").read_text())
    ref_dict = SeqIO.to_dict(SeqIO.parse("public_references/yeast-proteome.fasta", "fasta"))
    results = []
    for row in shortlist:
        orf = next(o for o in obj['unknown_orfs'] if o['location'] == row['location'])
        for candidate in row['candidates']:
            ref = ref_dict[candidate['id']]
            evidence = protein_evidence(orf['protein'], str(ref.seq))
            if evidence:
                results.append(dict(location=row['location'], reference_id=ref.id, description=ref.description, **evidence))
    with open("protein-evidence.json", "w") as f: json.dump(results, f, indent=2)
    print(json.dumps(results))

if __name__ == "__main__": main()