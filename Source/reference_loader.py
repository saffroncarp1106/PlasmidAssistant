from pathlib import Path
import os
from matching import load_database
from library_builder import build_library

def load_references(database_path):
    db_path = Path(database_path)
    if db_path.exists():
        return load_database(db_path)
    if db_path.name != "local_features.json":
        raise FileNotFoundError(db_path)

    snap_home = os.environ.get('SNAPGENE_HOME')
    roots = [db_path.parent.parent, db_path.parent.parent.parent]
    if snap_home:
        roots.append(Path(snap_home))

    entries = {}
    for root in roots:
        for name in ('Plasmids', 'sampleData'):
            folder = root / name
            if folder.is_dir():
                found, errors = build_library(folder)
                for entry in found:
                    entries.setdefault(entry['feature_id'], entry)

    return list(entries.values()) if entries else []
