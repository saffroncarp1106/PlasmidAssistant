# Plasmid annotation companion — development prototype

Open a sequence, find missing features, review matches, and save an annotated
copy. Built-in references load automatically. This uses local sequence matching
and seeded local alignment; it does not query BLAST or offer BLAST's coverage.

## Use

Keep the entire portable PlasmidAssistant folder together. Launch its executable,
choose Open sequence, then Find missing features. Review evidence and use Toggle
selected for suggestions you accept. Choose Save annotated copy and open that
new file in SnapGene. Existing annotations are preserved; existing files are
never overwritten. Delete the companion folder to remove it. Saved sequence
copies remain. The companion does not patch SnapGene or register an extension.

For .dna inputs, the window offers .dna output, preserving non-feature packets
byte for byte. GenBank export is also available but cannot preserve arbitrary
SnapGene metadata. FASTA/GenBank inputs export to GenBank. Set circular topology
for FASTA when appropriate.

Exact nonconflicting matches are selected by default. Variants, partial matches,
indels, conflicts and unknown ORFs require review. Evidence bands are heuristic,
not calibrated probabilities. Unknown ORFs do not establish gene identity or
expression. Only 37 references across 22 names are currently included. Absence
of a match does not mean a biological feature is absent.

## Development

Source requires Python and pip install -r requirements.txt.
Launch python review_ui.py, or:

    python annotate.py input.dna --output annotated.dna
    python annotate.py input.fasta --circular --output annotated.gb
    python -m unittest discover -s . -p "test_*.py"

CLI flags include --fast-only, --minimum-identity, --alignment-identity,
--alignment-coverage, --minimum-orf, --include-unknown-orfs, --accept FEATURE_ID
and --reject FEATURE_ID. Rejection overrides acceptance. A custom database is
optional, not required.

Local Ollama workers stage code for review during development. Packaged users
need neither Python nor Ollama. Historical dist, portable and release-native
builds differ; consult PROJECT_STATUS.md for current verification and next steps.

## Validation limits

Native tests preserve sequence, existing features and other packets. SnapGene
Viewer opened one generated 34,864 bp file with 12 original features and one
synthetic check feature. Broader documents, save/reload and clean Windows testing
remain pending. Analysis currently blocks the UI during searches. References
extracted from local files are not independently curated; redistribution rights
remain unresolved. See PROJECT_STATUS.md for benchmarks and completion checklist.
