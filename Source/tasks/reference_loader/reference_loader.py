from pathlib import Path
from matching import load_database
from library_builder import build_library

def load_references(database_path: Path):
    if database_path.exists(): return load_database(database_path)
    if database_path.name != 'local_features.json': raise FileNotFoundError(database_path)
    roots = [database_path.parent, database_path.parent.parent, database_path.parent.parent.parent]
    for root in roots:
        if not (root / 'Plasmids').exists() or not (root / 'sampleData').exists(): continue
        entries = []
        try: entries.extend(build_library(root))
        except Exception: pass
        if entries: return [e for e in entries if e['feature_id'] not in {x['feature_id'] for x in entries[:len(entries)//2]}]
    raise ValueError('No local reference examples were found. Place the companion folder inside the SnapGene installation containing Plasmids/sampleData, or choose a reference database under Advanced.')