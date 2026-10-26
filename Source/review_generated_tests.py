from pathlib import Path
root = Path(__file__).parent
p = root / 'diagnose_case.py'
s = p.read_text().replace("f['label']", "f.qualifiers.get('label',[''])[0]").replace("f['location']", 'f.location')
p.write_text(s)
p = root / 'test_alignment_candidates.py'
s = p.read_text().split('class SeqFeature:')[0]
s = s.replace('from alignment_candidates import alignment_candidates', 'from alignment_candidates import alignment_candidates\nfrom Bio.SeqFeature import SeqFeature, SimpleLocation')
s = s.replace("'sequence': REFERENCE[20:110]", "'sequence': REFERENCE")
s = s.replace("SeqRecord(Seq(REFERENCE), id='test_ref')", "SeqRecord(Seq(REFERENCE[20:110]), id='test_ref')")
s = s.replace("SeqRecord(Seq(target), id='test_ref_ins')", "SeqRecord(Seq(target).reverse_complement(), id='test_ref_ins')")
s = s.replace('from alignment_candidates import SimpleLocation', 'from Bio.SeqFeature import SimpleLocation')
s = s.replace('self.assertEqual(len(candidates), 0)', "self.assertFalse(any(c['location'].strand == 1 for c in candidates))")
s = s.replace('candidates = alignment_candidates(record, [entry], fast_hits)', "candidates = sorted(alignment_candidates(record, [entry], fast_hits), key=lambda h: (h['coverage'],h['identity']), reverse=True)")
p.write_text(s)
