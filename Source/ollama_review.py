"""Visible, logged local reviewer exchange; no automatic code execution."""
import json
import urllib.request
from pathlib import Path

root = Path(__file__).resolve().parent
prompt = '''Review this plasmid annotation implementation for concrete correctness bugs.
Focus on circular coordinates, strand handling, preserving existing annotations,
duplicate or conflicting matches, and safe export. Do not invent biological
identities or reference sequences. Return a short prioritized bug list with
reproduction examples, assumptions, and recommended fixes. You are a reviewer;
the primary agent makes scientific and architectural decisions.
'''
for name in ('annotate.py', 'matching.py'):
    prompt += '\nFILE: ' + name + '\n' + (root / name).read_text()
print('PRIMARY AGENT -> QWEN3.5:9B\n' + prompt, flush=True)
request = urllib.request.Request('http://localhost:11434/api/chat',
    data=json.dumps({'model': 'qwen3.5:9b', 'messages': [{'role': 'user', 'content': prompt}],
                     'stream': True, 'think': False}).encode(),
    headers={'Content-Type': 'application/json'})
print('\nQWEN3.5:9B -> PRIMARY AGENT\n', flush=True)
with (root / 'ollama-review.log').open('w', encoding='utf-8') as log:
    log.write(prompt + '\nRESPONSE:\n')
    try:
        with urllib.request.urlopen(request, timeout=600) as response:
            for line in response:
                chunk = json.loads(line)
                content = chunk.get('message', {}).get('content', '')
                print(content, end='', flush=True)
                log.write(content)
                log.flush()
    except Exception as error:
        print('\nReview failed:', error, flush=True)
        log.write('\nERROR: ' + str(error))
print('\n\nExchange complete. Log: ' + str(root / 'ollama-review.log'))
