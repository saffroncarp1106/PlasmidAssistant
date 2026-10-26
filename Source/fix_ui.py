from pathlib import Path
p = Path(__file__).with_name('review_ui.py')
text = p.read_text()
text = text.replace('self.chk_circular = ttk.Checkbutton', 'self.circular = tk.BooleanVar()\n        self.chk_circular = ttk.Checkbutton').replace('variable=tk.IntVar()', 'variable=self.circular')
text = text.replace('self.tree.bind("<<TreeviewSelect>>", self.on_select)', 'ttk.Button(control_frame, text="Toggle selected", command=self.on_select).grid(row=1, column=0)')
text = text[:text.index('    def open_sequence(self):')] + '''    def open_sequence(self):
        path = filedialog.askopenfilename(filetypes=[("Sequence", "*.fa *.fasta *.fna *.gb *.gbk *.genbank")])
        if not path:
            return
        try:
            self.record = read_record(path)
            self.input_path = Path(path)
            self.circular.set(self.record.annotations.get('topology') == 'circular')
            self.clear()
            self.root.title("Sequence Review: " + self.input_path.name)
        except Exception as error:
            messagebox.showerror("Open failed", str(error))

    def clear(self):
        self.row_to_feature.clear()
        self.selected_ids.clear()
        self.tree.delete(*self.tree.get_children())

    def load_database(self):
        path = filedialog.askopenfilename(filetypes=[("JSON", "*.json")])
        if path:
            try:
                load_database(Path(path))
                self.database_path = Path(path)
                self.clear()
            except Exception as error:
                messagebox.showerror("Database error", str(error))

    def toggle_topology(self):
        if self.record is not None and self.input_path.suffix.lower() in ('.fa', '.fasta', '.fna'):
            self.record.annotations['topology'] = 'circular' if self.circular.get() else 'linear'
            self.clear()
        elif self.record is not None:
            self.circular.set(self.record.annotations.get('topology') == 'circular')

    def add_row(self, feature, kind, evidence, selected=False, existing=False):
        name = feature.qualifiers.get('label', feature.qualifiers.get('gene', [feature.type]))[0]
        row = self.tree.insert('', tk.END, values=('Yes' if selected else '', name, kind, evidence, str(feature.location)))
        if not existing:
            self.row_to_feature[row] = feature
        if selected:
            self.selected_ids.add(row)

    def analyze(self):
        if self.record is None:
            return
        try:
            self.clear()
            _, candidates = annotate_sequence(self.record)
            entries = load_database(self.database_path) if self.database_path.exists() else []
            hits = match_features(self.record, entries)
            for feature in self.record.features:
                self.add_row(feature, 'Existing', 'Preserved', existing=True)
            for hit in hits:
                if hit['status'] in ('existing', 'superseded'):
                    continue
                selected = hit['confidence'] == 'High confidence' and hit['status'] == 'new'
                evidence = f"{hit['status']}; {hit['confidence']}; identity {hit['identity']:.1%}; coverage {hit['coverage']:.1%}"
                self.add_row(to_feature(hit), 'Known candidate', evidence, selected)
            for candidate in candidates:
                self.add_row(candidate, 'Unknown ORF', 'Identity unknown')
        except Exception as error:
            messagebox.showerror('Analysis failed', str(error))

    def on_select(self, event=None):
        for row in self.tree.selection():
            if row not in self.row_to_feature:
                continue
            if row in self.selected_ids:
                self.selected_ids.remove(row)
            else:
                self.selected_ids.add(row)
            self.tree.set(row, 'selected', 'Yes' if row in self.selected_ids else '')

    def export(self):
        if self.record is None:
            return
        path = filedialog.asksaveasfilename(defaultextension='.gb', initialfile=self.input_path.stem + '_autoannotated.gb')
        if not path:
            return
        try:
            import copy
            result = copy.deepcopy(self.record)
            result.features.extend(copy.deepcopy(self.row_to_feature[row]) for row in self.tree.get_children() if row in self.selected_ids)
            result.annotations.setdefault('molecule_type', 'DNA')
            stream = StringIO()
            SeqIO.write(result, stream, 'genbank')
            with open(path, 'x', encoding='utf-8') as handle:
                handle.write(stream.getvalue())
            messagebox.showinfo('Exported', path)
        except Exception as error:
            messagebox.showerror('Export failed', str(error))

if __name__ == '__main__':
    root = tk.Tk()
    App(root)
    root.mainloop()
'''
p.write_text(text)
