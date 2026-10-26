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
        orf_id = row.get("orf")
        if not orf_id: continue
        ref_seq = ref_dict.get(orf_id)
        if not ref_seq: continue
        evidence = protein_evidence(row["protein"], str(ref_seq.seq))
        results.append({**row, **evidence})
    with open("protein-evidence.json", "w") as f: json.dump(results, f, indent=2)
    print(json.dumps(results))

if __name__ == "__main__": main()