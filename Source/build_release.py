"""Rebuild a portable Windows folder without bundling development datasets."""
import shutil
import subprocess
import sys
from pathlib import Path

source = Path(__file__).resolve().parent
root = source.parent
subprocess.run([sys.executable, '-m', 'PyInstaller', '--noconfirm',
    '--distpath', str(root/'dist'), '--workpath', str(root/'.release-work'),
    str(source/'PlasmidAssistant.spec')], cwd=source, check=True)
target = root/'dist'/'PlasmidAssistant'
for name in ['protein_features.json','PROTEIN-REFERENCES.txt','README.md','START-HERE.txt']:
    shutil.copy2(root/name,target/name)
shutil.copytree(root/'licenses',target/'licenses',dirs_exist_ok=True)
print(target)
