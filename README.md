# PlasmidAssistant

A removable Windows companion for reviewing and adding plasmid annotations to
SnapGene `.dna` files. It runs in its own window and saves a new annotated copy.

## Run

1. Download and extract the repository ZIP, keeping `PlasmidAssistant.exe`,
   `_internal`, and `protein_features.json` together. Windows x64 is supported.
2. Run `PlasmidAssistant.exe`. Python and Ollama are not required.
3. Open a `.dna`, GenBank, or FASTA file and choose **Find missing features**.
4. Review candidates; double-click a row for evidence. Use **Toggle selected** to
   accept suggestions, then **Save annotated copy**. Open the new copy in SnapGene.

You can also drop a sequence file onto the executable in Windows Explorer.
Delete the extracted folder to uninstall. The program does not patch SnapGene,
register startup services, or overwrite the original sequence.

## References and BLAST

Bundled public protein references work without a SnapGene installation or manual
JSON loading. This is a small reference set, not a universal annotation database.
For additional DNA features, put the companion folder inside a SnapGene installation
containing `Plasmids` or `sampleData`; references are discovered automatically.
An optional `SNAPGENE_HOME` environment variable can point to that installation.
Without those DNA examples, the UI explicitly reports limited promoter/origin coverage.

**Search NCBI BLAST (online)** offers BLASTx against Swiss-Prot. It requests
confirmation before uploading the loaded DNA. It searches protein-coding regions,
not promoters. Internet access and NCBI availability are required; a search can
wait up to 15 minutes and can be cancelled. Cancel cannot retract an upload.
All BLAST candidates require review. Local search remains available offline.

Identity and reference coverage are shown separately; 70% identity is not a 70%
probability that a biological annotation is correct. Weak/partial matches and
unidentified ORFs are not automatically accepted. Protein searches use the standard
genetic code. Online queries are linear and may miss origin-spanning proteins;
local translated searches support circular and reverse-strand regions.

Promoter identification remains reference-dependent. This release does not claim
complete annotation of arbitrary plasmids or independent biological validation.

## Command line

```text
PlasmidAssistant.exe --cli input.dna --output annotated.dna
PlasmidAssistant.exe --cli input.fasta --circular --output annotated.gb
```

CLI output accepts exact, nonconflicting DNA matches by default. Use `--accept`
with a reviewed feature ID to include a weaker or protein candidate. Existing
output paths are refused. The GUI is recommended for inspecting candidates.

## Development

Tested with Python 3.14 on Windows. From the repository root:

```text
python -m pip install -r requirements-dev.txt
python -m unittest discover -s Source -p "test_*.py"
python Source/build_release.py
```

The build appears in `dist/PlasmidAssistant`; it contains no development datasets
or user plasmids. Source tests use synthetic fixtures. Historical scripts in
`Source` may refer to the original developer's local examples; use the commands
above for portable development and `Source/check_standalone.py` for distribution QA.
Public reference attribution is in `PROTEIN-REFERENCES.txt`; dependency licenses
are in `licenses` and `_internal`.

## Acknowledgment

Development included assistance from Qwen3.5 local agents and Codex, with generated
code reviewed and corrected before integration. The first archived release was
incomplete; this update fixes startup/search without installed DNA reference folders.
