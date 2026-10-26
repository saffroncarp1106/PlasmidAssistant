import argparse
import json
import time
import hashlib
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from matching import load_database
from benchmark_case import evaluate_record


def main():
    parser = argparse.ArgumentParser(description="Benchmark runner")
    parser.add_argument("corpus", type=Path)
    parser.add_argument("database", type=Path)
    parser.add_argument("--limit", type=int, default=20)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    entries = load_database(args.database)
    excluded_hashes = {e.get("source_sha256") for e in entries}

    groups = defaultdict(list)
    for path in sorted(args.corpus.rglob("*.dna")):
        groups[path.parent].append(path)

    ordered = []
    while any(groups.values()):
        for paths in groups.values():
            if paths:
                ordered.append(paths.pop(0))

    cases, errors = [], []
    excluded = 0

    for path in ordered:
        if len(cases) >= args.limit:
            break
        try:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            if digest in excluded_hashes:
                excluded += 1
                continue
            record = SeqIO.read(path, "snapgene")
            start = time.perf_counter()
            result = evaluate_record(record, entries)
            if result is not None:
                cases.append(
                    dict(
                        path=str(path),
                        bases=len(record),
                        seconds=time.perf_counter() - start,
                        **result
                    )
                )
        except Exception as e:
            errors.append(dict(path=str(path), error=str(e)))

    totals = {key: sum(c[key] for c in cases) for key in ["removed", "tp", "fp", "fn"]}
    if totals["tp"] + totals["fp"]:
        totals["precision"] = totals["tp"] / (totals["tp"] + totals["fp"])
    else:
        totals["precision"] = None
    if totals["tp"] + totals["fn"]:
        totals["recall"] = totals["tp"] / (totals["tp"] + totals["fn"])
    else:
        totals["recall"] = None

    report = dict(
        scope="Source-excluded annotation agreement, not biological validation",
        cases=cases,
        totals=totals,
        errors=errors,
        excluded_source_files=excluded
    )

    with args.output.open("x", encoding="utf-8") as f:
        f.write(json.dumps(report, separators=(",", ":")))
    print(json.dumps(totals, separators=(",", ":")))


if __name__ == "__main__":
    main()