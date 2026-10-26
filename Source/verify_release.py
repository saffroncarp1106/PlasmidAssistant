"""Execute release checks, then ask the local Ollama agent to audit real results."""
import json
import os
import subprocess
import sys
import urllib.request
from pathlib import Path

root=Path(__file__).resolve().parent.parent
checks=[['-m','unittest','discover','-s','Source','-p','test_*.py'],['Source/check_standalone.py']]
env=os.environ.copy()
env['PLASMID_DISTRIBUTION']=str(root/'dist'/'PlasmidAssistant')
results=[]
for args in checks:
    print('RUNNING',args,flush=True)
    try:
        run=subprocess.run([sys.executable,*args],cwd=root,env=env,capture_output=True,text=True,timeout=180)
        result={'command':args,'exit_code':run.returncode,'output':run.stdout+run.stderr}
    except subprocess.TimeoutExpired:
        result={'command':args,'exit_code':124,'output':'Timed out; not passed'}
    print(result['output'],flush=True)
    results.append(result)
(root/'release-verification.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
prompt='Audit these actual check results briefly. Say what passed, failed, and remains unproven. Never infer online BLAST success, GUI compatibility or biological accuracy from offline tests.\n'+json.dumps(results)
payload={'model':'qwen3.5:9b','messages':[{'role':'user','content':prompt}], 'think':False,'stream':True,
         'options':{'num_ctx':8192,'num_predict':500,'temperature':0.1}}
request=urllib.request.Request('http://localhost:11434/api/chat',data=json.dumps(payload).encode(),headers={'Content-Type':'application/json'})
print('OLLAMA AUDIT:',flush=True)
with (root/'release-verification-review.txt').open('w',encoding='utf-8') as log:
    with urllib.request.urlopen(request,timeout=600) as response:
        for line in response:
            text=json.loads(line).get('message',{}).get('content','')
            print(text,end='',flush=True)
            log.write(text)
            log.flush()
sys.exit(1 if any(r['exit_code'] for r in results) else 0)
