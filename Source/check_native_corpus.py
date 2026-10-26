from pathlib import Path
import tempfile, json
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation
from dna_writer import write_annotated_dna

root = Path(__file__).resolve().parent
sources = sorted((root.parent / 'Plasmids').rglob('*.dna'))[:30]
assert len(sources) == 30

with tempfile.TemporaryDirectory(prefix='native corpus ') as folder:
    for index, source in enumerate(sources):
        before = source.read_bytes()
        original = SeqIO.read(source, 'snapgene')
        assert len(original) >= 12
        synthetic = SeqFeature(SimpleLocation(2, 12, strand=-1), type='misc_feature', qualifiers={'label': ['Synthetic writer verification only']})
        output = Path(folder) / f'{index}.dna'
        write_annotated_dna(source, output, [synthetic])
        result = SeqIO.read(output, 'snapgene')
        assert result.seq == original.seq
        assert len(result.features) == len(original.features) + 1
        assert source.read_bytes() == before
        assert result.features[-1].location == synthetic.location

print(json.dumps({'verified_documents': len(sources), 'limitation': 'Native parser round trips only; no GUI or biological validation'}))