import tkinter as tk
import sys
from tkinter import ttk, filedialog, messagebox
from io import StringIO
from pathlib import Path
from Bio import SeqIO
from matching import load_database, match_features, to_feature
from reference_loader import load_references
from annotate import read_record, annotate_sequence

class App:
    def __init__(self, root):
        self.root = root
        self.root.title("Plasmid Auto-Annotation")
        self.record = None
        self.last_result = None
        self.busy = False
        self.cancel_job = None
        self.unknown_rows = []
        self.show_unknown = tk.BooleanVar(value=False)
        self.sensitive = tk.BooleanVar(value=True)
        self.status = tk.StringVar(value='Open a sequence to begin. References load automatically.')
        self.root.protocol('WM_DELETE_WINDOW', self.close)
        self.features = []  # List of (row_id, feature)
        self.selected_ids = set()
        self.row_to_feature = {}
        self.feature_to_row = {}
        self.database_path = (Path(sys.executable).parent if getattr(sys, 'frozen', False) else Path(__file__).parent) / "local_features.json"
        
        # UI Layout
        main_frame = ttk.Frame(root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Controls Frame
        control_frame = ttk.LabelFrame(main_frame, text="Controls", padding=5)
        control_frame.pack(fill=tk.X, pady=(0, 10))
        
        self.btn_open = ttk.Button(control_frame, text="Open sequence", command=self.open_sequence)
        self.btn_open.grid(row=0, column=0, padx=2, sticky='w')
        
        menu = tk.Menu(root)
        advanced = tk.Menu(menu, tearoff=False)
        advanced.add_command(label="Use custom reference database...", command=self.load_database)
        menu.add_cascade(label="Advanced", menu=advanced)
        root.config(menu=menu)
        ttk.Label(main_frame, text="References load automatically from your installation. Find missing features, review, then save a new annotated copy.").pack(fill=tk.X)
        
        self.btn_analyze = ttk.Button(control_frame, text="Find missing features", command=self.start_analysis)
        self.btn_analyze.grid(row=0, column=2, padx=2, sticky='w')
        
        self.circular = tk.BooleanVar()
        self.chk_circular = ttk.Checkbutton(control_frame, text="Circular topology (FASTA only)", variable=self.circular, command=self.toggle_topology)
        self.chk_circular.grid(row=0, column=3, padx=2, sticky='w')
        
        # Treeview Frame
        tree_frame = ttk.Frame(main_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True)
        
        columns = ("selected", "name", "kind", "evidence", "location")
        self.tree = ttk.Treeview(tree_frame, columns=columns, show="headings")
        self.tree.heading("selected", text="Selected")
        self.tree.heading("name", text="Name")
        self.tree.heading("kind", text="Kind")
        self.tree.heading("evidence", text="Evidence")
        self.tree.heading("location", text="Location (0-based)")
        
        self.tree.column("selected", width=80, anchor=tk.CENTER)
        self.tree.column("name", width=200)
        self.tree.column("kind", width=80)
        self.tree.column("evidence", width=150)
        self.tree.column("location", width=300)
        
        # Scrollbar
        scrollbar = ttk.Scrollbar(tree_frame, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)
        
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Export Button
        self.btn_export = ttk.Button(control_frame, text="Save annotated copy", command=self.export)
        self.btn_export.grid(row=0, column=4, padx=2, sticky='w')
        
        # Bind selection change
        ttk.Button(control_frame, text="Toggle selected", command=self.on_select).grid(row=1, column=0)
        ttk.Checkbutton(control_frame, text='Show unknown ORFs (not identified genes)', variable=self.show_unknown, command=self.filter_unknown).grid(row=1, column=2, columnspan=2)
        self.btn_cancel = ttk.Button(control_frame, text='Cancel search', command=self.cancel_search, state='disabled')
        self.btn_cancel.grid(row=1, column=4)
        self.chk_sensitive = ttk.Checkbutton(control_frame, text='Sensitive suggestions (70%+ identity; review required)', variable=self.sensitive, command=self.clear)
        self.chk_sensitive.grid(row=2, column=0, columnspan=5, sticky='w')
        self.btn_online = ttk.Button(control_frame, text='Search NCBI BLAST (online)', command=self.start_online)
        self.btn_online.grid(row=3, column=0, columnspan=3, sticky='w')
        self.tree.bind('<Double-1>', self.show_evidence)
        ttk.Label(main_frame, textvariable=self.status, wraplength=900).pack(fill=tk.X)
        
    def open_sequence(self, path=None):
        if self.busy:
            return
        path = path or filedialog.askopenfilename(filetypes=[("Sequence", "*.dna *.fa *.fasta *.fna *.gb *.gbk *.genbank")])
        if not path:
            return
        try:
            input_bytes = Path(path).read_bytes()
            record = read_record(path)
            if Path(path).read_bytes() != input_bytes:
                raise ValueError('Source changed while loading; reopen it.')
            self.record = record
            self.input_bytes = input_bytes
            self.input_path = Path(path)
            self.circular.set(self.record.annotations.get('topology') == 'circular')
            self.clear()
            self.root.title("Sequence Review: " + self.input_path.name)
            self.status.set(f'{self.input_path.name}: {len(record):,} bp. Ready to find missing features.')
        except Exception as error:
            messagebox.showerror("Open failed", str(error))

    def clear(self):
        self.last_result = None
        for row in self.unknown_rows:
            if self.tree.exists(row):
                self.tree.delete(row)
        self.unknown_rows.clear()
        self.row_to_feature.clear()
        self.selected_ids.clear()
        self.tree.delete(*self.tree.get_children())

    def load_database(self):
        if self.busy:
            return
        path = filedialog.askopenfilename(filetypes=[("JSON", "*.json")])
        if path:
            try:
                load_database(Path(path))
                self.database_path = Path(path)
                self.clear()
            except Exception as error:
                messagebox.showerror("Database error", str(error))

    def toggle_topology(self):
        if self.busy:
            return
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
        if kind == 'Unknown ORF':
            self.unknown_rows.append(row)

    def filter_unknown(self):
        for row in self.unknown_rows:
            if self.show_unknown.get():
                self.tree.move(row, '', 'end')
            else:
                self.selected_ids.discard(row)
                self.tree.set(row, 'selected', '')
                self.tree.detach(row)

    def set_busy(self, busy):
        self.busy = busy
        for widget in (self.btn_open, self.btn_analyze, self.btn_export, self.chk_circular, self.chk_sensitive, self.btn_online):
            widget.configure(state='disabled' if busy else 'normal')
        self.btn_cancel.configure(state='normal' if busy else 'disabled')

    def start_analysis(self):
        if self.record is None or self.busy:
            return
        from copy import deepcopy
        from background_job import start_job
        from analysis_service import analyze_record, AnnotationConfig
        snapshot, database = deepcopy(self.record), self.database_path
        config = AnnotationConfig(alignment_identity=.7 if self.sensitive.get() else .8)
        self.clear()
        self.set_busy(True)
        self.status.set('Searching references... You can cancel this search.')
        def work():
            return analyze_record(snapshot, load_references(database), config)
        def success(result):
            self.cancel_job = None
            self.set_busy(False)
            self.render_result(result)
        def failure(error):
            self.cancel_job = None
            self.set_busy(False)
            self.status.set('Search failed. Check the reference database and try again.')
            messagebox.showerror('Analysis failed', str(error))
        self.cancel_job = start_job(self.root, work, success, failure)

    def cancel_search(self):
        if self.cancel_job:
            self.cancel_job()
            self.cancel_job = None
        self.set_busy(False)
        self.status.set('Search cancelled. No annotations were saved.')

    def close(self):
        if self.cancel_job:
            self.cancel_job()
        self.root.destroy()

    def analyze(self):
        if self.record is None:
            return
        try:
            self.clear()
            entries = load_references(self.database_path)
            from analysis_service import analyze_record
            result = analyze_record(self.record, entries)
            self.render_result(result)
        except Exception as error:
            messagebox.showerror('Analysis failed', str(error))

    def render_result(self, result):
        self.clear()
        self.last_result = result
        candidates, hits = result['unknown_orfs'], result['known_features']
        for feature in self.record.features:
            self.add_row(feature, 'Existing', 'Preserved', existing=True)
        for hit in hits:
            if hit['status'] in ('existing', 'superseded'):
                continue
            selected = hit['confidence'] == 'High confidence' and hit['status'] == 'new'
            units = 'AA' if hit.get('method') in ('protein_local_alignment', 'blastx') else 'DNA'
            evidence = f"{hit['status']}; {hit['confidence']}; {units} identity {hit['identity']:.1%}; coverage {hit['coverage']:.1%}"
            self.add_row(to_feature(hit), 'Known candidate', evidence, selected)
        for candidate in candidates:
            self.add_row(candidate, 'Unknown ORF', 'Identity unknown')
        self.filter_unknown()
        self.status.set(f'Review ready: {len(self.row_to_feature) - len(candidates)} known candidates; {len(candidates)} unknown ORFs hidden unless enabled. Save creates a new copy.')

    def on_select(self, event=None):
        for row in self.tree.selection():
            if row not in self.row_to_feature:
                continue
            if row in self.selected_ids:
                self.selected_ids.remove(row)
            else:
                self.selected_ids.add(row)
            self.tree.set(row, 'selected', 'Yes' if row in self.selected_ids else '')

    from native_export import export
    from online_review import start_online

    def show_evidence(self, event=None):
        row = self.tree.focus()
        feature = self.row_to_feature.get(row)
        if feature is not None:
            messagebox.showinfo('Annotation evidence', '\n\n'.join(
                [self.tree.set(row, 'name'), str(feature.location)] + feature.qualifiers.get('note', [])))

if __name__ == '__main__':
    if '--cli' in sys.argv:
        sys.argv.remove('--cli')
        # Windowed executable has no console stdout.
        import io
        if sys.stdout is None:
            sys.stdout = io.StringIO()
        from annotate import main
        main()
    else:
        root = tk.Tk()
        app = App(root)
        if len(sys.argv) > 1:
            root.after(0, lambda: app.open_sequence(sys.argv[1]))
        root.mainloop()
