import json
import shutil
import tempfile
from pathlib import Path
from reference_loader import load_references

import unittest

class ReferenceLoaderTests(unittest.TestCase):
    def test_real_automatic_reference_discovery(self):
        # Read local_features.json relative to this file
        source_file = Path(__file__).parent / "local_features.json"
        
        with tempfile.TemporaryDirectory() as tmp_dir:
            install_path = Path(tmp_dir)
            
            # Create directory structure
            plasmids_dir = install_path / "Plasmids"
            companion_dir = install_path / "Companion"
            plasmids_dir.mkdir()
            companion_dir.mkdir()
            
            # Copy source file to Plasmids/example.dna (simulating the expected location)
            source = Path(json.loads(source_file.read_text())[0]['source'])
            shutil.copyfile(source, plasmids_dir / 'example.dna')
            
            # Call load_references with a non-existent local_features.json in Companion
            result = load_references(companion_dir / "local_features.json")
            
            # Assert result is non-empty
            assert len(result) > 0
            
            # Assert all entries have sequence and feature_id
            for entry in result:
                assert 'sequence' in entry
                assert 'feature_id' in entry
                
            # Assert reference IDs (feature_ids) are unique
            feature_ids = [entry['feature_id'] for entry in result]
            assert len(feature_ids) == len(set(feature_ids))
            
            # Assert no local_features.json was created in Companion
            assert not (companion_dir / "local_features.json").exists()
    
    def test_missing_custom_filename_raises_file_not_found(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            install_path = Path(tmp_dir)
            companion_dir = install_path / "Companion"
            companion_dir.mkdir()
            
            # Try to load a custom filename that doesn't exist and isn't 'local_features.json'
            try:
                load_references(companion_dir / "custom_name.json")
                assert False, "Should have raised FileNotFoundError"
            except FileNotFoundError:
                pass # Expected behavior
