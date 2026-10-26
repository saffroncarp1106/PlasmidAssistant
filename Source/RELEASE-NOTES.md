# Release notes — 2026-10-24

The removable local annotation companion is built in this folder. Launch
PlasmidAssistant.exe; keep _internal and the other files alongside it.

Implemented: native .dna/GenBank/FASTA input, preservation of existing features,
automatic references from installed examples, exact/variant/partial/indel matching,
circular and reverse-strand handling, unknown ORFs, background search/cancellation,
review and competing-selection guards, native .dna copies, file arguments, and
removal by deleting the folder. The program performs no online sequence search.

The local qwen3.5:9b development worker generated implementation and verification
code, with corrections before integration, and audited actual verification output.
45 automated tests passed. Thirty native documents round-tripped; the relocated
executable worked in a path with spaces/Unicode and empty PATH using automatically
discovered example references. A deliberately removed origin was recovered by the
actual executable. Original source bytes and annotations were preserved.

The companion GUI loaded and searched a real 34,864 bp native file. SnapGene
Viewer previously opened a generated native copy showing all 12 original features
plus one synthetic writer-check annotation. These observations do not establish
all SnapGene versions, normal save/reload, or a clean Windows installation.

Biological limitations: references come from your installed examples, not an
independently curated universal database. Searches use heuristic evidence bands
and bounded alignment windows, and can miss features. Unassigned ORFs are not gene
identifications. Source-excluded annotation agreement was imperfect (historical
reports remain in the development source). Do not treat this as validated BLAST
coverage or calibrated biological certainty.

Future validation: independently curated broader references; held-out biological
positive/negative controls; clean Windows and additional SnapGene versions. These
are not claimed complete by this release. No proprietary reference database is
included in this folder; the program reads the recipient's installed examples.
