import tempfile
import tkinter as tk
import unittest
from pathlib import Path
from unittest.mock import patch
from Bio import SeqIO
from main import App
from test_engine import record


class ReviewTests(unittest.TestCase):
    def test_review_export_and_overwrite_protection(self):
        root = tk.Tk()
        root.withdraw()
        self.addCleanup(root.destroy)
        app = App(root)
        app.record = record('ATG' + 'AAA' * 100 + 'TAA')
        app.input_path = Path('fixture.fa')
        app.database_path = Path(__file__).with_name('local_features.json')
        app.show_unknown.set(True)
        app.analyze()
        rows = app.tree.get_children()
        self.assertTrue(rows)
        self.assertFalse(app.selected_ids)
        app.tree.selection_set(rows[0])
        app.on_select()
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / 'result.gb'
            with patch('main.filedialog.asksaveasfilename', return_value=str(target)), patch('main.messagebox.showinfo'), patch('main.messagebox.showerror') as error:
                app.export()
                error.assert_not_called()
                output = SeqIO.read(target, 'genbank')
                self.assertEqual(len(output.features), 1)
                self.assertEqual(len(app.record.features), 0)
                before = target.read_bytes()
                app.export()
                error.assert_called_once()
                self.assertEqual(before, target.read_bytes())


if __name__ == '__main__':
    unittest.main()
