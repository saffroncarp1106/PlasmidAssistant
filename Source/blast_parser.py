import sys
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.SeqFeature import SimpleLocation

def blast_candidates(xml_text, record, rid):
    """
    Parses BLASTX XML results and filters candidates based on specific criteria.
    
    Args:
        xml_text (str): XML string from NCBI BLAST result.
        record (Bio.SeqRecord): The query sequence record.
        rid (str): RID parameter for the source URL.
        
    Returns:
        list: List of dictionaries containing filtered hit information.
    """
    
    parsed = NCBIXML.read(StringIO(xml_text))
    if parsed.application.upper() != "BLASTX" or parsed.query_length != len(record):
        raise ValueError("BLAST program or query length mismatch")

    results = []

    for alignment in parsed.alignments:
        if parsed.application.upper() != 'BLASTX':
            raise ValueError("Application must be BLASTX")
        
        if parsed.query_length != len(record.seq):
            raise ValueError("Query length mismatch")

        for hsp in alignment.hsps:
            align_len = hsp.align_length
            if align_len <= 0 or alignment.length <= 0:
                continue
                
            identity = hsp.identities / align_len
            subj_start = hsp.sbjct_start
            subj_end = hsp.sbjct_end
            
            length_on_subject = abs(subj_end - subj_start) + 1
            
            coverage = length_on_subject / alignment.length

            if not (0.7 <= identity <= 1 and 0.7 <= coverage <= 1 and 0 <= hsp.expect <= 1e-5):
                continue
            
            
            try:
                q_start = hsp.query_start
                q_end = hsp.query_end
            except AttributeError:
                continue

            min_q = min(q_start, q_end)
            max_q = max(q_start, q_end)
            
            if min_q < 1 or max_q > len(record.seq):
                continue
                
            try:
                frame = hsp.frame[0]
                if frame not in (-3,-2,-1,1,2,3):
                    continue
                if frame > 0:
                    strand = 1
                else:
                    strand = -1
            except (IndexError, TypeError):
                continue

            loc = SimpleLocation(min_q - 1, max_q, strand=strand)

            entry = {
                'feature_id': f'blast:{alignment.accession}',
                'canonical_name': alignment.hit_def,
                'feature_type': 'misc_feature',
                'source': f'https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID={rid}',
                'version': rid,
            }
            hit = {
                'entry': entry,
                'location': loc,
                'identity': identity,
                'coverage': coverage,
                'match': 'protein_variant',
                'confidence': 'Review required',
                'status': 'new',
                'method': 'blastx',
                'evalue': hsp.expect,
                'bitscore': hsp.bits
            }
            
            results.append(hit)

    return results
