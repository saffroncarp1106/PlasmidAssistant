"""Explicitly requested online protein identification, outside the Tk event loop."""
def start_online(self):
    from copy import deepcopy
    from threading import Event
    from tkinter import messagebox
    from background_job import start_job
    from blast_transport import BlastClient
    from blast_parser import blast_candidates
    from analysis_service import analyze_record
    from reference_loader import load_references
    from resolver import resolve_hits

    if self.record is None or self.busy:
        return
    if not messagebox.askyesno('Search NCBI BLAST?',
            'Send this DNA sequence to NCBI for BLASTx against Swiss-Prot? '
            'This identifies protein-coding regions; it does not identify promoters. '
            'Results require review. The search may take several minutes.'):
        return
    snapshot = deepcopy(self.record)
    prior = deepcopy(getattr(self, 'last_result', None))
    database = self.database_path
    cancelled = Event()
    self.set_busy(True)
    self.status.set('NCBI BLASTx: searching Swiss-Prot. Cancel stops waiting and future requests.')

    def work():
        rid, xml = BlastClient(cancel=cancelled).search(str(snapshot.seq))
        result = prior or analyze_record(snapshot, load_references(database))
        result['known_features'] = resolve_hits(result['known_features'] + blast_candidates(xml, snapshot, rid))
        return result

    def success(result):
        self.cancel_job = None
        self.set_busy(False)
        self.render_result(result)
        self.status.set('BLAST search complete. Review candidates and their evidence before saving a copy.')

    def failure(error):
        self.cancel_job = None
        self.set_busy(False)
        self.status.set('BLAST search failed; existing review results are preserved.')
        messagebox.showerror('NCBI BLAST', str(error))

    discard = start_job(self.root, work, success, failure)
    def cancel():
        cancelled.set()
        discard()
    self.cancel_job = cancel
