from pathlib import Path
from matching import load_database
from library_builder import build_library


def load_references(database_path):
    database_path = Path(database_path)
    if database_path.exists():
        return load_database(database_path)
    if database_path.name != 'local_features.json':
        raise FileNotFoundError(database_path)
    for root in (database_path.parent.parent, database_path.parent.parent.parent):
        entries = {}
        for name in ('Plasmids', 'sampleData'):
            folder = root / name
            if folder.is_dir():
                found, errors = build_library(folder)
                for entry in found:
                    entries.setdefault(entry['feature_id'], entry)
        if entries:
            return list(entries.values())
    raise ValueError('No local reference examples were found. Place the companion folder inside the SnapGene installation containing Plasmids/sampleData, or choose a reference database under Advanced.')
