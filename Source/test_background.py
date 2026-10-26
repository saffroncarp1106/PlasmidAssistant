import unittest
import tkinter as tk
import threading
import time
from background_job import start_job

class TestBackgroundJob(unittest.TestCase):
    def setUp(self):
        self.root = tk.Tk()
        self.root.withdraw()
        self.results = []
        self.main_id = threading.get_ident()
        self.gate = threading.Event()
        self.work = lambda: (self.gate.wait(2), 42)[1]

    def tearDown(self):
        self.root.update()
        self.root.destroy()

    def test_callback_after_gate(self):
        cancel = start_job(self.root, self.work, lambda value: self.results.append((value, threading.get_ident())), self.fail)
        self.assertEqual(len(self.results), 0)
        self.gate.set()
        deadline = time.monotonic() + 3
        while not self.results and time.monotonic() < deadline:
            self.root.update()
            time.sleep(0.01)
        self.root.update()
        self.assertEqual(self.results, [(42, self.main_id)])

    def test_cancel_immediately(self):
        cancel = start_job(self.root, self.work, lambda value: self.results.append((value, threading.get_ident())), self.fail)
        cancel()
        self.gate.set()
        time.sleep(0.1)
        self.root.update()
        self.assertEqual(len(self.results), 0)

if __name__ == '__main__':
    import tkinter as tk
    unittest.main()