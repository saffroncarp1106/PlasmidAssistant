# Portable SnapGene companion requirement

The deliverable is a self-contained folder that can be placed inside a SnapGene
installation or any writable location. A bundled executable launches a separate
review window. Users remove the companion by deleting its folder.

Do not patch SnapGene executables, DLLs, registry entries, startup configuration,
or bundled data. Do not require administrator privileges. Runtime dependencies
must be bundled for end users. Ollama is a development worker, not an end-user
dependency or a source of biological identities.

Current integration is file based: select a sequence, review results, save a native
.dna copy for .dna inputs (or GenBank), then open/reload it in SnapGene. File drop
is still pending. Native output preserves non-feature packets. Never write GenBank
bytes under a .dna extension. Any future replacement workflow must be explicit, create
a backup, detect changes to the source since loading, and save atomically.

Packaging acceptance: launch on a clean Windows machine without Python; run from
a path containing spaces; operate without network; preserve the installation and
inputs; uninstall by deletion; verify dependency and reference-data redistribution
licenses before bundling. These checks remain pending.
