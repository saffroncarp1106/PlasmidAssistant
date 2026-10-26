import threading
from queue import Empty, Queue
import tkinter as tk

def start_job(root, work, on_success, on_error):
    q = Queue()
    cancelled = [False]
    
    def worker():
        try:
            result = work()
            q.put((True, result))
        except Exception as e:
            q.put((False, e))
    
    t = threading.Thread(target=worker, daemon=True)
    t.start()
    
    after_id_list = [None]
    
    def poll():
        after_id_list[0] = None
        if cancelled[0]: return
        try:
            success, data = q.get_nowait()
        except Empty:
            try:
                after_id_list[0] = root.after(50, poll)
            except tk.TclError:
                cancelled[0] = True
            return
        if success: on_success(data)
        else: on_error(data)
    
    def cancel():
        cancelled[0] = True
        if after_id_list[0] is not None:
            try:
                root.after_cancel(after_id_list[0])
            except tk.TclError:
                pass
        after_id_list[0] = None
    
    after_id_list[0] = root.after(0, poll)
    return cancel
