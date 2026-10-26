"""Run fixed project checks and ask the local worker to audit their actual output."""
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

root = Path(__file__).resolve().parent
checks = [['-m', 'unittest', 'discover', '-s', '.', '-p', 'test_*.py'], ['check_example_protein.py'], ['check_objective.py']]
results = []
for check in checks:
    print('RUNNING:', ' '.join(check), flush=True)
    try:
        run = subprocess.run([sys.executable, *check], cwd=root, capture_output=True, text=True, timeout=180)
        item = dict(command=check, exit_code=run.returncode, output=run.stdout + run.stderr)
    except subprocess.TimeoutExpired:
        item = dict(command=check, exit_code=124, output='Timed out after 180 seconds; NOT passed')
    results.append(item)
    print(item['output'], flush=True)
(root / 'optimization-verification-results.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
prompt = ('Audit these actual verification results. Report what passed, what failed, and what they do NOT prove. '
          'Do not claim SnapGene GUI compatibility or independent accuracy from unit tests. Be concise.\n' + json.dumps(results))
request = urllib.request.Request('http://localhost:11434/api/chat', data=json.dumps({
    'model': 'qwen3.5:9b', 'messages': [{'role': 'user', 'content': prompt}],
    'think': False, 'stream': True, 'options': {'num_ctx': 8192, 'num_predict': 600, 'temperature': 0.1}}).encode(),
    headers={'Content-Type': 'application/json'})
print('QWEN VERIFICATION REVIEW:', flush=True)
with (root / 'optimization-verification-review.txt').open('w', encoding='utf-8') as log:
    with urllib.request.urlopen(request, timeout=1200) as response:
        for line in response:
            content = json.loads(line).get('message', {}).get('content', '')
            print(content, end='', flush=True)
            log.write(content)
            log.flush()
