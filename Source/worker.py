"""Visible local implementation worker; writes only a named staging file."""
import argparse
import json
import urllib.request
from pathlib import Path

root = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('task')
args = parser.parse_args()
task = json.loads((root / 'tasks' / (args.task + '.json')).read_text(encoding='utf-8'))
prompt = task['prompt']
for filename in task.get('context', []):
    prompt += '\nFILE ' + filename + '\n' + (root / filename).read_text(encoding='utf-8')
print('PRIMARY -> QWEN\n' + prompt + '\nQWEN -> PRIMARY', flush=True)
schema = {'type': 'object', 'properties': {'code': {'type': 'string'}, 'report': {'type': 'string'}}, 'required': ['code', 'report'], 'additionalProperties': False}
request = urllib.request.Request('http://localhost:11434/api/chat', data=json.dumps({
    'model': 'qwen3.5:9b', 'messages': [{'role': 'user', 'content': prompt}],
    'think': task.get('think', False), 'stream': True,
    'options': {'temperature': 0.1, 'num_ctx': 8192, 'num_predict': task.get('num_predict', 4096)}}).encode(), headers={'Content-Type': 'application/json'})
stage = root / 'tasks' / args.task
stage.mkdir(exist_ok=True)
chunks = []
with (stage / 'exchange.log').open('w', encoding='utf-8') as log:
    log.write(prompt + '\nRESPONSE\n')
    with urllib.request.urlopen(request, timeout=1200) as response:
        for line in response:
            message = json.loads(line).get('message', {})
            thinking = message.get('thinking', '')
            if thinking:
                print(thinking, end='', flush=True)
                log.write(thinking)
                log.flush()
            content = message.get('content', '')
            chunks.append(content)
            print(content, end='', flush=True)
            log.write(content)
            log.flush()
raw = ''.join(chunks).strip()
if raw.startswith('```'):
    raw = raw.split('\n', 1)[1].rsplit('```', 1)[0].strip()
try:
    result = json.loads(raw)
except json.JSONDecodeError:
    result = {'code': raw, 'report': 'Raw source response'}
code = result.get('code', result.get('corrected_code'))
if not isinstance(code, str) or not code.strip():
    raise ValueError('Worker did not return nonempty source code')
compile(code, task['output'], 'exec')
with (stage / task['output']).open('x', encoding='utf-8') as handle:
    handle.write(code)
print('\nSTAGED FOR REVIEW\n' + result['report'], flush=True)
