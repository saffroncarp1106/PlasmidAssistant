# Project status and handoff — 2026-10-24

## Current deliverable

`../PlasmidAssistant/PlasmidAssistant.exe` is the current removable companion.
The older dist, portable and release-candidate directories are historical builds.
Keep the entire new folder together. START-HERE.txt and RELEASE-NOTES.md in it
explain usage, requirements, and validation limits.

Completed local functionality:
- Native .dna, GenBank and FASTA input; existing annotation preservation.
- Reference discovery automatically from the recipient's installed Plasmids and
  sampleData. No extracted reference database is redistributed in the new folder.
  Development local_features.json remains a local test library.
- Exact and full-length variant matching, seeded partial/indel suggestions,
  reverse/circular handling, complete unknown ORFs with linear scan performance.
- Background GUI search, cancellation, progress/status, unknown-ORF filtering,
  file arguments (Explorer drop onto executable), explicit evidence review.
- Equal-score boundary ambiguities require review. Competing selected features
  cannot be saved together. Native output preserves non-feature packets, refuses
  overwrite and changed-source exports. GenBank export remains available.
- Self-contained executable runtime, dependency notices, no registry changes,
  installation patching, background service, persistent reference cache or network.

Verification is run through visible local Ollama development workers. Generated
code is reviewed before execution because qwen3.5:9b frequently invents APIs or
weakens assertions. Actual outputs: verification-results.json; local model audit:
verification-review.txt. Current suite: 45 automated tests. Additional scripts
cover 30 native documents, native overwrite refusal, relocated delivery with an
empty PATH and automatic references, and recovery of a deliberately removed ori
using the actual executable. All passed at the final build checkpoint.

GUI inspection: built companion opened and searched native-ui-check.dna, showed
search status and preserved existing annotations while hiding unknown ORFs.
SnapGene Viewer opened the earlier generated file with 12 original features plus
one synthetic writer-check annotation. Do not generalize this to untested versions.

## Still required for a broadly validated biological product

- Independently curate reference identities, provenance, aliases and appropriate
  boundaries; broaden yeast auxotrophic markers, origins and regulatory elements.
  Installed sample coverage varies. The delivery requires installed examples or
  an explicitly chosen custom reference JSON; it does not contain a universal DB.
- Expand held-out biological benchmarks and calibrate thresholds. Existing
  benchmark-source-excluded.json: 30 plasmids, 93/109 recovered, 18 unmatched,
  16 missed using the old 26-reference library. benchmark-library37.json before
  the boundary fix: 91/109, 20 unmatched, 18 missed on the same cases. Do not
  mislabel source annotation agreement as independent biological truth.
- Evaluate false identities separately from boundary differences and missing
  original annotations. Synthetic negatives do not prove biological specificity.
- Verify a genuinely clean Windows machine, protected install locations, more
  SnapGene versions and native save/reload with primers/history/unknown packets.
- Alignment windows are bounded and can miss divergent/repeated features. Protein
  comparison and explicit opt-in BLAST escalation are future extensions.
- UI can be improved with richer grouping, configurable GUI thresholds and a
  dedicated drop target. Current CLI already exposes thresholds; Explorer can
  drop files onto the executable. Cancel discards results, not computation.

## Rebuild and verify

Source: annotation-assistant. `python -m PyInstaller --noconfirm --distpath ..
--workpath build-delivery PlasmidAssistant.spec` rebuilds the owned companion
folder. Close the GUI first. Reapply START-HERE/RELEASE-NOTES and dependency license
files after rebuilding; PyInstaller replaces its output directory. Do NOT copy
local_features.json into the redistributable folder. Keep source/installation
plasmids untouched. `python -u verify_local.py` runs real checks and has Ollama
audit their output; show this in a visible cmd window. Task specs and transcripts
are in tasks/. Rejected model outputs are not shipped.
