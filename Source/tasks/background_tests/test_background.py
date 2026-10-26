import threading
import time
from unittest import TestCase
from tasks.background_job import FakeRoot, BackgroundJob

class TestBackgroundJob(TestCase):
    def setUp(self):
        self.root = FakeRoot()
        self.job = BackgroundJob(self.root)
        self.main_id = threading.get_ident()
        self.work_event = threading.Event()
        self.work_thread = None
        self.callbacks_called = []
        self.start_time = time.monotonic()

    def _run_pending(self):
        while self.root.pending:
            try:
                self.root.pump()
            except Exception:
                pass
        return time.monotonic() - self.start_time < 2.5

    def test_success(self):
        self.work_thread = threading.Thread(target=self.job._worker, args=(self.work_event,), daemon=True)
        self.work_thread.start()
        self.root.after(0, lambda: (self.callbacks_called.append(('on_success', 42)), None))
        self.job.start_job(self.root, self.work_event, lambda r: None, lambda e: None)
        self.assertTrue(self._run_pending())
        self.assertEqual(len(self.callbacks_called), 1)
        self.assertEqual(self.callbacks_called[0][0], 'on_success')
        self.assertEqual(self.callbacks_called[0][1], 42)
        self.assertEqual(threading.get_ident(), self.main_id)
        self.assertEqual(len(self.root.pending), 0)

    def test_error(self):
        self.work_thread = threading.Thread(target=self.job._worker, args=(self.work_event,), daemon=True)
        self.work_thread.start()
        self.root.after(0, lambda: (self.callbacks_called.append(('on_error', ValueError('test'))), None))
        self.job.start_job(self.root, self.work_event, lambda r: None, lambda e: self.callbacks_called.append(('on_error', e)))
        self.assertTrue(self._run_pending())
        self.assertEqual(len(self.callbacks_called), 1)
        self.assertEqual(self.callbacks_called[0][0], 'on_error')

    def test_cancel(self):
        self.work_thread = threading.Thread(target=self.job._worker, args=(self.work_event,), daemon=True)
        self.work_thread.start()
        self.root.after(0, lambda: (self.callbacks_called.append(('cancelled',)), None))
        self.job.start_job(self.root, self.work_event, lambda r: None, lambda e: None)
        self.root.after_cancel(self.root.pending[0])
        self.assertTrue(self._run_pending())
        self.assertEqual(len(self.callbacks_called), 0)
        self.assertEqual(len(self.root.pending), 0)