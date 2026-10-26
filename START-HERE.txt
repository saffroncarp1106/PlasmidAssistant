Plasmid Annotation Companion - local-first release

1. Keep this entire PlasmidAssistant folder together inside your SnapGene install.
   The install must contain its Plasmids and/or sampleData examples.
2. Run PlasmidAssistant.exe. No Python or Ollama installation is required.
3. Open a .dna, GenBank or FASTA sequence, then click Find missing features.
   You may also drop a sequence file onto the executable in Windows Explorer.
4. Review candidates. Exact, nonconflicting matches are selected automatically.
   Toggle selected accepts/rejects other candidates. Competing selections cannot
   be saved together. Unknown ORFs are hidden and unselected by default.
5. Save annotated copy. For .dna input this creates a new native .dna file.
   Open that copy in SnapGene. Original files are not overwritten.
6. Optional: Search NCBI BLAST (online) sends the loaded DNA to NCBI only after
   confirmation. BLASTx searches Swiss-Prot for protein-coding regions. It does not
   identify promoters or guarantee every gene is represented. Internet is required.
   Double-click a candidate to inspect identity, reference coverage and provenance.
7. To uninstall, close the program and delete this folder. No registry changes,
   SnapGene patches, startup hooks or external cache are created.

References are built automatically in memory from YOUR installed examples.
This distribution does not include extracted SnapGene reference sequences.
Public protein references in protein_features.json supplement the installed DNA
examples. Their source URLs and attribution are included in PROTEIN-REFERENCES.txt.
No samples? The program explains where to place the folder; Advanced also allows
an explicit custom reference JSON. Coverage depends on the installed examples.

Search runs in the background. Cancel discards the result; the worker may finish
in the background, but never writes annotations. Closing exits the companion.
FASTA topology must be set correctly. Location notation is zero-based, half-open.
Evidence bands are heuristic, not calibrated biological probabilities. Partial,
variant and conflicting matches require review; unknown ORFs are not gene IDs.
Local search uses pairwise sequence comparison. The optional online action uses
actual NCBI BLASTx. BLAST results are suggestions, not automatic biological truth.
Online searches use the standard genetic code and linear query coordinates;
origin-spanning proteins may be missed online. Local searches handle circular DNA.
Cancel stops future online polls; it cannot retract a sequence already submitted.

Verification and known limits are in RELEASE-NOTES.md. Native format preservation
and local functionality were tested; clean-Windows and independent biological
accuracy validation are not established. Keep normal research backups.

Command line:
  PlasmidAssistant.exe --cli input.dna --output annotated.dna
  PlasmidAssistant.exe --cli input.fasta --circular --output annotated.gb

Dependencies: Python 3.14.7, Biopython 1.88, NumPy 2.5.4, Tcl/Tk, and the
PyInstaller 6.22.3 bootloader. License notices are in licenses and the NumPy
metadata licenses under _internal. No affiliation with SnapGene is claimed.
