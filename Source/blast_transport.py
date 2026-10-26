import time
import threading
import urllib.parse
import urllib.request
import re
from typing import Tuple, Optional

class BlastClient:
    def __init__(self, cancel: Optional[threading.Event] = None, timeout: int = 900, last_request: float = 0.0):
        self.cancel = cancel or threading.Event()
        self.timeout = timeout
        self.last_request = last_request

    def _request(self, params: dict) -> str:
        
        # Check cancellation with a small sleep window to avoid busy-waiting too aggressively
        wait_time = max(0.0, 10.0 - (time.monotonic() - self.last_request))
        if wait_time > 0:
            self.cancel.wait(wait_time)
            
        if self.cancel.is_set():
            raise RuntimeError("Request cancelled")

        url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
        data = urllib.parse.urlencode(params).encode('utf-8')
        
        req = urllib.request.Request(url, data=data, method='POST')
        # Set a reasonable timeout for the network call itself
        self.last_request = time.monotonic()
        with urllib.request.urlopen(req, timeout=30) as response:
            return response.read().decode('utf-8')

    def search(self, sequence: str) -> Tuple[str, str]:
        sequence = str(sequence).upper()
        # Validate input: nonempty IUPAC DNA, max 100k bases
        if not sequence:
            raise ValueError("Sequence cannot be empty")
        
        # Basic validation for IUPAC DNA characters (A, C, G, T, R, Y, S, W, K, M, B, D, H, V, N)
        valid_chars = set("ACGTRYSWKMBDHVN")
        if not all(c in valid_chars for c in sequence):
            raise ValueError("Sequence contains invalid IUPAC DNA characters")
        
        if len(sequence) > 100000:
            raise ValueError("Sequence exceeds maximum length of 100,000 bases")

        # Prepare submission parameters
        params = {
            "CMD": "Put",
            "PROGRAM": "blastx",
            "DATABASE": "swissprot",
            "QUERY": sequence,
            "EXPECT": "1e-5",
            "HITLIST_SIZE": "50",
            "TOOL": "PlasmidAssistant"
        }

        # Submit request
        response_text = self._request(params)
        
        # Extract RID using regex
        rid_pattern = r"RID\s*=\s*([A-Z0-9-]+)"
        match = re.search(rid_pattern, response_text)
        
        if not match:
            raise RuntimeError("No RID found in response")
            
        rid = match.group(1).strip()

        # Polling loop for results
        deadline = time.monotonic() + self.timeout
        poll_interval = 60.0
        
        while True:
            current_time = time.monotonic()
            
            # Check cancellation
            if self.cancel.is_set():
                raise RuntimeError("Polling cancelled")
                
            # Check deadline
            if current_time >= deadline:
                raise TimeoutError(f"Search timed out after {self.timeout} seconds")

            # Wait before next poll (cancellation-aware)
            wait_time = min(poll_interval, deadline - current_time)
            self.cancel.wait(wait_time)
            
            if self.cancel.is_set():
                raise RuntimeError("Polling cancelled")

            if time.monotonic() >= deadline:
                raise TimeoutError("BLAST search deadline reached")

            # Fetch status
            status_params = {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_TYPE": "XML"
            }
            
            try:
                poll_response = self._request(status_params)
            except RuntimeError as e:
                if str(e) == "Request cancelled":
                    raise
                # Re-raise other runtime errors immediately
                raise

            # Check for completion or failure
            if "<BlastOutput" in poll_response:
                return (rid, poll_response)
                
            if "Status" in poll_response and ("FAILED" in poll_response or "UNKNOWN" in poll_response):
                # Handle optional whitespace around = e.g., Status = FAILED
                status_match = re.search(r'Status\s*=\s*(\w+)', poll_response)
                if status_match:
                    status_val = status_match.group(1).strip()
                    if status_val in ("FAILED", "UNKNOWN"):
                        raise RuntimeError(f"BLAST search failed with status: {status_val}")

        # Should not be reached due to deadline check or return/raise above
        raise TimeoutError("Search did not complete within timeout")