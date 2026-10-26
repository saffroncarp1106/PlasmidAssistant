import unittest
from pathlib import Path
from tempfile import TemporaryDirectory, NamedTemporaryFile
from unittest.mock import Mock, patch
from types import SimpleNamespace
from native_export import export


class TestExportGuard(unittest.TestCase):
    def test_source_change_guard(self):
        with TemporaryDirectory() as tmp:
            source_path = Path(tmp) / "source.dna"
            source_path.write_bytes(b'old')
            
            obj = SimpleNamespace(
                record=Mock(),
                input_path=source_path,
                input_bytes=b'old',
                selected_ids=set(),
                row_to_feature={},
                tree=Mock()
            )
            obj.tree.get_children.return_value = []
            
            with patch('tkinter.filedialog.asksaveasfilename', return_value=str(Path(tmp) / "out.dna")), \
                 patch('tkinter.messagebox.showerror') as error_mock, \
                 patch('tkinter.messagebox.showinfo'), \
                 patch('native_export.dna_writer.write_annotated_dna') as writer_mock:
                
                export(obj)
                
                error_mock.assert_called_once()
                self.assertIn(b'changed', error_mock.call_args.args[1])
                writer_mock.assert_not_called()
                self.assertFalse((Path(tmp) / "out.dna").exists())
                self.assertEqual(source_path.read_bytes(), b'changed')


if __name__ == '__main__':
    unittest.main()