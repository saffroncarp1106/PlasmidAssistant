# Project handoff and completion checklist

Updated 2026-10-24. This is an unfinished prototype, not a validated annotation product.

## User's required delivery

A removable, self-contained folder placed inside a SnapGene installation. Launch
a companion window, open a .dna file, automatically find missing features using
references loaded automatically, review evidence, save an annotated .dna copy,
then reopen/reload in SnapGene. No SnapGene binary patches, no manual database
setup, no Python/Ollama requirement for the end user, no overwritten originals.

Development must predominantly use local Ollama qwen3.5:9b. Show exchanges in
visible cmd windows. Keep primary-model tokens low. Review generated code before
integration: this model has repeatedly invented APIs, ignored constraints, and
produced invalid tests. Passing its self-review is insufficient evidence.

## Current implementation

- Source lives here, separate from the surrounding SnapGene binary installation.
- GenBank/FASTA/.dna reading, preserved existing annotations, circular ORFs.
- Exhaustive seeded full-length matching on both strands; 450 comparisons against
  exhaustive scoring passed. ~0.10 seconds for 26 references on a 34,864 bp plasmid.
- Seeded local alignment for partial/indel suggestions; review-only; finite search
  windows can miss divergent/repeated features. Integrated into analysis_service,
  CLI and UI source. ~0.8 seconds with 37 references on the same plasmid.
- Current local_features.json contains 37 sequence references across 22 names,
  extracted from installation Plasmids and sampleData. Provenance hashes recorded.
  Local source annotations are not independently curated truth. Redistribution
  rights for this reference data have not been established.
- Native writer preserves raw non-feature packets, appends feature XML, reparses
  before saving, refuses overwrite. Tested simple/reverse/circular features and a
  real file with 12 existing annotations. Installed SnapGene Viewer opened
  native-ui-check.dna and displayed all 13 features (12 original plus a clearly
  labelled synthetic writer check). Save/reload and broader GUI checks remain.
- Review UI auto-loads references and supports annotated .dna copies. GenBank
  remains available. CLI also supports native output. Source-change guards exist.
- Linear ORF implementation integrated. 1,800 deterministic oracle comparisons
  cover linear/circular scans, ambiguity and sequence lengths modulo three.
- Latest source suite: 35 tests passed, including boundary ambiguity and changed
  source export refusal tests. Ollama generated the new verification programs;
  primary review corrected faulty assertions/imports before integration.
  Equal-score same-name matches with different boundaries now require review
  instead of selecting an arbitrary database ID. A failed source load also retains
  the previously loaded record rather than partially replacing its state.
- Current build: release-candidate/PlasmidAssistant/PlasmidAssistant.exe, with
  local_features.json beside it. Rebuilt after ORF and boundary fixes. Older dist,
  portable and release-native folders are historical. Frozen verification results
  are recorded by verify_local.py; consult verification-results.json.
- Current candidate passed native executable smoke verification, preserving the
  sequence and all 12 original features/qualifiers. The explicit unknown-ORF test
  added 222 unassigned ORFs (234 total); these are NOT identified genes.
- Ollama's native corpus verifier passed 30 installation documents. This is parser
  round-trip verification, not SnapGene GUI verification or biological validation.
- The candidate also passed after copying its complete folder to a path containing
  spaces and Greek alpha, with PATH empty and PYTHONHOME/PYTHONPATH removed. This
  confirms relocation on this machine, not operation on a clean Windows machine.
- Ollama audited the actual outputs: verification-review.txt. All five commands
  in verification-results.json returned exit code 0. No verification task remains
  running at this checkpoint. release-candidate/build-hashes.json records exe and
  reference-library SHA256 values.

## Benchmark evidence

benchmark-source-excluded.json: 30 plasmids, 109 removed features, 93 recovered,
18 unmatched predictions, 16 misses. Strict agreement precision .8378, recall
.8532. This excludes reference-source file hashes but is NOT independent biological
validation. Five unmatched calls differ in origin boundaries; 13 are promoter
matches absent from supplied annotations. Do not relabel these as true positives
without independent evidence. benchmark-discrepancies.json records the audit.
These results used the earlier 26-reference library; repeat with current source
and record database/code versions. Do not overwrite ground truth or prior reports.

benchmark-library37.json (before boundary fix): same 30 plasmids, 109 removed,
91 recovered, 20 unmatched, 18 missed. Two regressions occur in pET SUMO
(linearized) and pETBlue-1 AccepTor. This exposed arbitrary same-name boundary
tie-breaking. The fix now retains these ambiguities for review; rerun and audit
the benchmark rather than tuning it to match stored annotations. This benchmark
exercises fast matching, not the entire partial/indel pipeline.

## Remaining before calling the project done

Next priority: responsive background analysis and review UX, then independent
reference curation/negative benchmarks. Current GUI searches still block the
window. Preserve the existing candidate while implementing the next version.

- [x] Integrate and verify linear ORF span implementation against an independent
  exhaustive oracle, including n modulo 3, ambiguity, reverse and origin crossing.
- [ ] Regression-test analysis_service end to end, configuration validation,
  partial/indel exports, source-change refusal, duplicates, conflicts, and negatives.
- [ ] Resolve candidate deduplication and competing reference boundaries. Avoid
  declaring high-confidence functional identity solely from a short sequence motif.
  Keep evidence bands configurable and explicitly uncalibrated.
- [ ] Broaden and independently curate reference sequences, provenance, aliases,
  feature classes, and lengths. Yeast HIS3/LEU2/TRP1/ADE2/MET15 and many origins,
  promoters and terminators remain absent from the small local library.
- [ ] Build a held-out benchmark spanning feature classes and exact, mutated,
  truncated, reversed, circular, overlapping, and unknown constructs. Report
  identity/class errors separately from boundary errors and missing source labels.
  Establish acceptance criteria and evaluate specificity on negative controls.
- [ ] Validate native output in the actual SnapGene UI across multiple versions/
  documents. Confirm features, colors, primer/history/unknown packet preservation,
  reload behavior and normal subsequent saves. Biopython read-back alone is not this.
- [ ] Keep UI responsive with background analysis, meaningful status, actionable
  errors, understandable coordinates, grouped existing/new/unknown results, and
  selected-feature conflict handling. Unknown ORFs currently clutter large results.
- [ ] Finish open/drop-file or executable argument workflow. File-based integration
  is implemented; no supported desktop extension hook has been established.
- [ ] Rebuild latest source with dependencies and library, package notices/licenses,
  remove development artifacts from the shipped folder, and update user docs.
- [ ] Verify frozen end-to-end native annotation, clean Windows without Python,
  offline operation, paths with spaces/Unicode, protected install locations,
  and removal by deleting the companion folder without changing SnapGene.
- [ ] Update DELIVERY and finalize release instructions. README and candidate
  START-HERE now describe the actual native export and automatic reference workflow.
- [ ] Optional later: translated protein comparison and external BLAST escalation.
  Do not upload user sequences without an explicit search choice.

## How to continue locally

Task specs: tasks/*.json. Runner: `python -u worker.py TASK_NAME` in visible cmd.
Outputs and logs are staged under tasks/TASK_NAME; never auto-trust them.
`python -u verify_local.py` runs fixed checks and asks Ollama to audit actual results.
`python -m unittest discover -s . -p "test_*.py"` runs source tests here.
`python -m PyInstaller --noconfirm --distpath RELEASE_FOLDER PlasmidAssistant.spec`
builds a folder executable. Copy the CURRENT local_features.json beside its exe.
Keep original installation files and source plasmids unchanged.
