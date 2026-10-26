from dataclasses import dataclass
from copy import deepcopy
import matching
import resolver
import alignment_candidates

@dataclass
class AnnotationConfig:
    minimum_orf: int = 300
    minimum_identity: float = .95
    search_partial: bool = True
    alignment_identity: float = .7
    alignment_coverage: float = .5
    max_windows: int = 20
    search_protein: bool = True

    def validate(self):
        if self.minimum_orf < 6:
            raise ValueError("minimum_orf must be at least 6")
        if self.max_windows < 1:
            raise ValueError("max_windows must be at least 1")
        for attr in ["minimum_identity", "alignment_identity", "alignment_coverage"]:
            val = getattr(self, attr)
            if not (0 < val <= 1):
                raise ValueError(f"{attr} must be in (0, 1]")

def analyze_record(record, entries, config=None):
    if config is None:
        config = AnnotationConfig()
    config.validate()
    
    from annotate import annotate_sequence
    preserved, unknown = annotate_sequence(record, config.minimum_orf)
    hits = matching.match_features(record, entries, config.minimum_identity)
    
    if config.search_partial:
        extended_hits = alignment_candidates.alignment_candidates(
            record, entries, hits,
            config.alignment_identity, config.alignment_coverage, config.max_windows
        )
        hits.extend(extended_hits)
    
    if config.search_protein:
        import json
        import sys
        from pathlib import Path
        from protein_candidates import protein_candidates
        base = Path(sys.executable).parent if getattr(sys, 'frozen', False) else Path(__file__).parent
        reference_file = base / 'protein_features.json'
        if reference_file.exists():
            references = json.loads(reference_file.read_text(encoding='utf-8'))
            hits.extend(protein_candidates(record, unknown, references))
    resolved_hits = resolver.resolve_hits(hits)
    
    unknown = [orf for orf in unknown if not any(
        hit['status'] == 'new' and hit['confidence'] == 'High confidence'
        and hit['location'].strand == orf.location.strand
        and set(orf.location).issubset(set(hit['location'])) for hit in resolved_hits)]
    
    return {
        "record": preserved,
        "dna_reference_count": len(entries),
        "existing_features": deepcopy(record.features),
        "known_features": resolved_hits,
        "unknown_orfs": unknown
    }
