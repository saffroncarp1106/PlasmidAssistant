import threading
from types import SimpleNamespace
import unittest
from unittest.mock import patch, MagicMock, call
from io import StringIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from blast_parser import blast_candidates
from blast_transport import BlastClient

class TestBlastCandidates(unittest.TestCase):
    
    def setUp(self):
        self.mock_record = SeqRecord(Seq('A'*300), id='test_query')
        self.valid_xml = 'mock parsed separately'
        self.hsp = SimpleNamespace(align_length=100,identities=90,sbjct_start=1,sbjct_end=100,expect=1e-20,bits=80,query_start=300,query_end=1,frame=(-1,0))
        alignment = SimpleNamespace(length=100,accession='X',hit_def='Example',hsps=[self.hsp])
        self.parsed = SimpleNamespace(application='BLASTX',query_length=300,alignments=[alignment])
        patcher = patch('blast_parser.NCBIXML.read',return_value=self.parsed)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_valid_hit_parsing(self):
        """Test parsing a valid BLASTX hit with correct strand and identity."""
        result = blast_candidates(self.valid_xml, self.mock_record, "test_rid")
        
        # Assert exactly one hit is returned
        self.assertEqual(len(result), 1)
        
        hit = result[0]
        
        # Assert location details: strand -1, start 0 (1-based to 0-based conversion logic in code), end 300
        self.assertEqual(hit['location'].strand, -1)
        self.assertEqual(hit['location'].start, 0)
        self.assertEqual(hit['location'].end, 300)
        
        # Assert identity is approx 90% (allowing float precision)
        self.assertAlmostEqual(hit['identity'], 0.9, places=2)
        
        # Assert entry details
        self.assertEqual(hit['entry']['feature_id'], 'blast:X')
        self.assertEqual(hit['entry']['canonical_name'], 'Example')
        self.assertEqual(hit['method'], 'blastx')

    def test_invalid_application_raises(self):
        """Test that non-BLASTX application raises ValueError."""
        self.parsed.application = "BLASTN"
        xml_bad = self.valid_xml
        with self.assertRaises(ValueError) as context:
            blast_candidates(xml_bad, self.mock_record, "test_rid")
        self.assertIn("BLAST program or query length mismatch", str(context.exception))

    def test_wrong_query_length_raises(self):
        """Test that mismatched query length raises ValueError."""
        short_record = SeqRecord(Seq('A'*100), id='short')
        with self.assertRaises(ValueError) as context:
            blast_candidates(self.valid_xml, short_record, "test_rid")
        self.assertIn("BLAST program or query length mismatch", str(context.exception))

    def test_invalid_xml_raises(self):
        """Test that invalid XML raises ValueError via mock side_effect."""
        bad_xml = "<invalid>xml</invalid>"
        
        with patch('blast_parser.NCBIXML.read') as mock_read:
            mock_read.side_effect = ValueError("Invalid XML structure")
            
            with self.assertRaises(ValueError):
                blast_candidates(bad_xml, self.mock_record, "test_rid")

    def test_weak_identity_filtered(self):
        """Test that hits with low identity (e.g., 60%) are filtered out."""
        # Modify XML to have 60% identity
        self.hsp.identities = 60
        xml_low_id = self.valid_xml
        
        result = blast_candidates(xml_low_id, self.mock_record, "test_rid")
        self.assertEqual(len(result), 0)

    def test_transport_invalid_dna_rejected(self):
        """Test that BlastClient rejects invalid DNA before making network calls."""
        client = BlastClient()
        
        # Test with invalid character 'X'
        with patch.object(client, '_request') as mock_request:
            mock_request.return_value = "RID=TEST123"
            
            with self.assertRaises(ValueError):
                client.search("ACGTXYZ")
            
            # Ensure network was not called
            mock_request.assert_not_called()

    def test_transport_cancel_before_urlopen(self):
        """Test that a pre-set Event cancels the request before urlopen is called."""
        cancel_event = threading.Event()
        cancel_event.set() # Pre-set
        
        client = BlastClient(cancel=cancel_event)
        
        with patch('urllib.request.urlopen') as mock_urlopen:
            with self.assertRaises(RuntimeError):
                client.search("ACGT")
            
            mock_urlopen.assert_not_called()

if __name__ == '__main__':
    import threading
    unittest.main()