def export(self):
    from copy import deepcopy
    from pathlib import Path
    import io
    import dna_writer
    from tkinter import filedialog, messagebox
    if self.record is None:
        return
    native = self.input_path.suffix.lower() == '.dna'
    default_ext = '.dna' if native else '.gb'
    filetypes = [('SnapGene DNA', '*.dna'), ('GenBank', '*.gb')] if native else [('GenBank', '*.gb')]
    path = filedialog.asksaveasfilename(
        defaultextension=default_ext,
        initialfile=self.input_path.stem + '_autoannotated' + default_ext,
        filetypes=filetypes
    )
    if not path:
        return
    try:
        selected = [deepcopy(self.row_to_feature[row]) for row in self.tree.get_children() if row in self.selected_ids]
        from selection_guard import validate_selection
        validate_selection(selected)
        suffix = Path(path).suffix.lower()
        if suffix == '.dna':
            if not native:
                raise ValueError('Native format required for .dna files.')
            if self.input_path.read_bytes() != self.input_bytes:
                raise ValueError('Source file changed since loading.')
            dna_writer.write_annotated_dna(self.input_path, Path(path), selected)
        elif suffix in ('.gb', '.gbk', '.genbank'):
            from Bio import SeqIO
            record_copy = deepcopy(self.record)
            record_copy.features.extend(selected)
            record_copy.annotations['molecule_type'] = 'DNA'
            stream = io.StringIO()
            SeqIO.write(record_copy, stream, 'genbank')
            with Path(path).open('x', encoding='utf-8') as f:
                f.write(stream.getvalue())
        else:
            raise ValueError(f'Unsupported extension: {suffix}')
        messagebox.showinfo('Exported', f'New annotated copy saved to:\n{path}')
    except Exception as error:
        messagebox.showerror('Export failed', str(error))
