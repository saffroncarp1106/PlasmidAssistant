import threading
from queue import Queue
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
    
    def poll():
        if cancelled[0]: return
        try:
            success, data = q.get_nowait()
            if success: on_success(data)
            else: on_error(data)
        except Exception: pass
        
        try:
            root.after(50, poll)
        except tk.TclError:
            pass
    
    def cancel():
        cancelled[0] = True
        try: root.after_cancel(poll) if hasattr(root.after, 'cancel') else None
        except tk.TclError: pass
    
    poll()
    return cancel