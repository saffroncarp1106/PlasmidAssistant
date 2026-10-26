import random
import unittest
from orf_spans import complete_orf_spans


def naive(seq, circular):
    scan = seq * 2 if circular else seq
    n = len(seq)
    output = []
    for start in range(n):
        if scan[start:start+3] != 'ATG':
            continue
        limit = start + n if circular else n
        for end in range(start+3, limit-2, 3):
            codon = scan[end:end+3]
            if set(codon) - set('ACGT'):
                break
            if codon in ('TAA', 'TAG', 'TGA'):
                output.append((start, end+3-start))
                break
    return sorted(output)


class TestOrfSpans(unittest.TestCase):
    
    def test_naive_reference(self):
        seq = 'ATGATGAAATAA'
        expected = [(0, 12), (3, 9)]
        self.assertEqual(naive(seq, False), expected)
        self.assertEqual(naive(seq, True), expected)

    def test_random_linear(self):
        lengths = [9, 10, 11, 30, 31, 32, 60, 61, 100]
        for length in lengths:
            seq = ''.join(random.choices('ACGTN', k=length))
            result = complete_orf_spans(seq, False, 6)
            expected = naive(seq, False)
            self.assertEqual(result, expected)

    def test_random_circular(self):
        lengths = [9, 10, 11, 30, 31, 32, 60, 61, 100]
        for length in lengths:
            seq = ''.join(random.choices('ACGTN', k=length))
            result = complete_orf_spans(seq, True, 6)
            expected = naive(seq, True)
            self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()