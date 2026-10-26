# Optimization in progress (2026-10-24)

User requests minimal primary tokens, predominantly local Ollama implementation
AND verification, visible cmd windows. Original release remains under
../PlasmidAssistant with its dated ZIP. Current source has UNRELEASED changes.

User supplied C:/Users/llkin/Documents/objective and amount of plasmid annotation.pdf.
Read/rendered its one page: assignment expects 7 features, including revealed
largest ORF; promoters/terminators/antibiotic+auxotrophic markers. It gives NO names.
Source target: C:/Users/llkin/Documents/WhoamI plasmid new_DAB.dna, 6574bp,
existing ori [634:1223](-) and CYC1 terminator [4780:5028](+).
New example: C:/Users/llkin/Documents/plasmid example.dna, 6550bp DIFFERENT sequence.
user-example-features.json holds annotations. AmpR, URA3 (reverse in example!),
mRFP1 [4105:4780](+) starts without ATG, ori, NOS terminator, CYC1 terminator,
URA3 promoter and a primer. Never copy example coordinates to original blindly.

Source improvements so far:
- library_builder expands installed-label extraction from 37/22 to 128 refs/82 names
  (expanded-library-optimization.json). Promoters/terminators selected by labels;
  provenance remains installation annotations, not independent curation.
- seed_windows now ranks seed-supported windows; .7 identity sensitive mode uses
  8mer seeds step4. analysis defaults .7 DNA identity/.5 coverage for suggestions.
- GUI sensitive checkbox, explicit AA vs DNA identity; never probability claims.
- alignment_candidates skips references already exact-found (may miss extra partial
  copies), prunes overlapping same-reference candidates by coverage, and suppresses
  partial origin suggestions wholly inside already annotated origins.
- protein_evidence.py reuses corrected local alignment coordinate algorithm.
- protein_candidates.py matches translated ATG ORFs to protein_features.json;
  all protein hits review-only. analysis_service integrates optional file besideexe.
- NEXT: translated_segments.py worker pending, removes ATG requirement so mRFP1
  boundary can be recovered (stop-to-stop six-frame protein regions). Integrate
  careful strand/circular coordinate mapping and test before rebuilding.

Evidence (no user sequence was uploaded anywhere):
- Downloaded UniProt public yeast proteome, then reviewed proteins length300..700
  (68MB gzip) to public_references. Local 4mer shortlist plus protein alignment.
- Largest ORF [2755:4105](+) matches Candida albicans ASH1 UniProt Q5A5Q6:
  449aa, 99.5546% identity, 100% reference/query coverage with standard genetic code.
  Standard vs organism-specific genetic-code differences have NOT been evaluated.
- [5508:6321](+) is URA3-like: 71.863% AA identity, 97.75% reference coverage to
  S. cerevisiae P03962. Do not claim it is specifically S. cerevisiae URA3.
- mRFP1 CAH64892.1 public NCBI reference: current ATG-only search [4135:4780](+),
  100% AA identity, 95.11% reference coverage. Example starts at4105; fix segment search.
- Exact DNA new AmpR [1394:2255](-), NOS terminator [6321:6574](+).
- Original-promoter identification still unresolved. Do NOT force total7 or claim
  full assignment completion merely because 2 existing+5 new feature names=7.
- User answered yes to BLAST permission, then provided example and said goal is
  autoannotation from BLAST. No remote BLAST run yet; current matcher isn't BLAST.
  Need explicit BLAST provider/import with provenance if fulfilling literal scope.

Files: objective-final-report.json, objective-analysis.json, protein-shortlist*.json,
protein-evidence.json, user-example-features.json, protein_features.json (8 public
references). UniProt-license.json proves CC BY4.0; preserve attribution for bundle.
Sources: https://www.uniprot.org/uniprotkb/Q5A5Q6/entry ;
https://www.ncbi.nlm.nih.gov/protein/CAH64892.1 ;
https://www.fpbase.org/protein/mrfp1/ ; UniProt license https://www.uniprot.org/help/license

Ollama generated verification scripts needed substantial corrections. Reject
self-assertions and invented APIs. inspect_objective.py, assess_proteins.py,
check_objective.py corrected; check_objective asserts three protein identities
and unchanged source bytes, writes objective-final-report.json. Earlier verify_local
output still45tests/stale executable because new tests weren't all run/integrated.
Need run updated suite; update check_known_recovery's exact-feature-count assertion
since expanded library can legitimately add features to original annotated input.
Check preservation by prefix/qualifiers instead. Negative controls must use expanded
library and protein refs. Rebuild named optimization release only after checks.
Copy protein_features.json beside exe and include public-source/license attribution.
Keep inputs unchanged, export only new clearly labelled candidate copies.

Some long combined exec commands were rejected with only 'blocked by policy'.
Reviewed, smaller read-only downloads and explicit file patches succeeded. No
user data was sent to any external service. Local Ollama prompts contain excerpts
and evidence as authorized. Do not circumvent any approval denial.


## 2026-10-25 continuation
Ollama was offline; restarted visibly (`ollama serve`). Local model staged
segments_fix, protein_segments, segment_tests, example_check. Drafts required
correction: missed stop advancement/size guard, invented location dictionaries,
and vacuous assertions. Segment test draft was rejected (reimplemented target).
Integrated translated_segments.py, corrected protein_candidates.py six-frame
stop-to-stop search, and actual imported-function tests. check_example_protein.py
passed: example mRFP1 [4105:4780](+), identity 1.0, coverage .9955556, source bytes
unchanged. check_objective.py passes with ASH1, mRFP1, URA3 evidence. verify_local.py
now includes both actual-objective checks and sends real outputs to Ollama.
Review verification-results.json timestamps and exit codes before release.
Current frozen executable is STILL the old build; source optimization unreleased.
Remaining: stronger reverse/circular/fusion coordinate tests; expanded-library
negative controls; promoter identification; real opt-in BLAST integration;
protein-reference attribution and packaging; new frozen build and GUI/export QA.
Do not claim project complete or local pairwise alignment is BLAST.

Verified continuation: optimization-verification-results.json has 52 passing unit tests plus passing example and original-objective checks. Ollama completed optimization-verification-review.txt. Its phrase '10 distinct feature matches' counts candidate reports, not 10 validated biological features. Expanded full-suite negative-control runtime remains under investigation. New executable not built.
