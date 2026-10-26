# Validation — 2026-10-26

- 65 source tests passed on Windows with Python 3.14.7.
- Rebuilt executable ran from a temporary folder containing only the executable,
  runtime and bundled protein references, with an empty PATH and no installed
  SnapGene examples. A public-reference protein annotation was exported to native
  DNA format; sequence and existing annotation were preserved.
- The executable processed the developer's supplied real plasmid with installed
  DNA examples enabled. It preserved its two existing features and exported five
  additional reviewed/exact matches. The source bytes stayed unchanged. This is
  a regression check, not independent evidence of complete biological annotation.
- The packaged GUI launched and displayed its controls. Full end-to-end GUI use
  across SnapGene versions and a clean Windows installation are not established.
- The local Ollama agent audited actual unit/standalone check output. Generated
  implementation drafts were corrected before integration.

The executable hash is recorded in build-manifest.json. Online BLAST depends on
NCBI availability; offline checks do not establish successful online searches.
Promoter coverage remains reference-dependent. No claim of universal or complete
plasmid annotation is made.
