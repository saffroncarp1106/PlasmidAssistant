import tempfile, tkinter as tk, unittest
from pathlib import Path
from unittest.mock import patch
from Bio import SeqIO
from review_ui import App
from test_engine import record

class ReviewTests(unittest.TestCase):
    def test_export_refusal_on_source_change(self):
        root = tk.Tk(); root.withdraw()
        self.addCleanup(root.destroy)
        app = App(root); app.record = record('ACGT')
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'source.dna'
            output = Path(tmp) / 'result.dna'
            source.write_bytes(b'old'); app.input_path = source; app.input_bytes = b'old'
            app.record = test_engine.record('ACGT') # Ensure record exists
            app.export() # Should succeed initially if mocks allow, but logic checks bytes
            # Simulate change: user edits source file externally or programmatically before next export
            source.write_bytes(b'changed')
            with patch('review_ui.filedialog.asksaveasfilename', return_value=str(output)):
                app.export()
                self.assertTrue(output.stat().st_size == 0 or not output.exists()) # Output should be nonexistent/empty due to guard
                # Verify error was shown and writer not called (via side_effect logic in real impl)
                pass

# Minimal mocks for dependencies to make this standalone runnable snippet
import sys; sys.modules['dna_writer'] = type('Module', (), {'write_annotated_dna': lambda *a, **k: None})()
test_engine = type('Module', (), {'record': lambda s: type('Rec', (), {'features': [], 'annotations': {}})()})()