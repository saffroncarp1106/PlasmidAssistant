from Bio.Seq import Seq

def translated_segments(sequence, circular=False, minimum_aa=80):
    sequence = str(sequence).upper()
    n = len(sequence)
    results = []
    stops = {'TAA', 'TAG', 'TGA'}
    
    for strand, oriented in [(1, sequence), (-1, str(Seq(sequence).reverse_complement()))]:
        scan = oriented * 2 if circular else oriented
        for frame in range(3):
            start = frame
            pos = frame
            while pos < len(scan) - 2:
                codon = scan[pos:pos+3]
                if codon in stops or set(codon) - set('ACGT'):
                    if pos - start >= minimum_aa * 3 and start < n:
                        results.append({
                            'start': start, 
                            'strand': strand, 
                            'protein': str(Seq(scan[start:pos]).translate()),
                            'has_stop': codon in stops
                        })
                pos += 3
            
            end = frame + ((len(scan) - frame) // 3) * 3
            if end - start >= minimum_aa * 3 and start < n:
                results.append({
                    'start': start, 
                    'strand': strand, 
                    'protein': str(Seq(scan[start:end]).translate()),
                    'has_stop': False
                })
    return results