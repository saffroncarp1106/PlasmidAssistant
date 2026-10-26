import sys, json, pathlib, Bio.SeqIO, collections
d = pathlib.Path(sys.argv[1])
c = collections.Counter()
e = []
for p in d.rglob('*.dna'):
    try:
        r = Bio.SeqIO.read(p, 'snapgene')
        c.update(f.qualifiers.get('label',[''])[0] for f in r.features)
    except Exception as x:
        e.append({'file': str(p), 'error': str(x)})
print(json.dumps({'labels': c.most_common(), 'errors': e, 'files': len(list(d.rglob('*.dna')))}))
