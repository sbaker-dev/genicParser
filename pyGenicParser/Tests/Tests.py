from .. import *

from pathlib import Path
import unittest


class MyTestCase(unittest.TestCase):

    @staticmethod
    def test_data_path():
        """Validate all our paths"""
        assert Path(Path(__file__).parent, "Data", "EUR.ldpred_21.bed").exists()
        assert Path(Path(__file__).parent, "Data", "EUR.ldpred_21.bim").exists()
        assert Path(Path(__file__).parent, "Data", "EUR.ldpred_21.fam").exists()
        assert Path(Path(__file__).parent, "Data", "EUR.ldpred_21.bgen").exists()
        assert Path(Path(__file__).parent, "Data", "EUR.ldpred_21.bgen.bgi").exists()
        assert Path(Path(__file__).parent, "Data", "Write", "write.txt").exists()

    @staticmethod
    def test_bgen_bgi_write():
        """Test writing bgen.bgi"""
        bgen_path = Path(Path(__file__).parent, "Data", "EUR.ldpred_21.bgen")
        write_path = Path(Path(__file__).parent, "Data", "Write")
        BgenObject(bgen_path, bgi_write_path=write_path).create_bgi()

        out_path = Path(Path(__file__).parent, "Data", "Write", "EUR.ldpred_21.bgen.bgi")
        assert out_path.exists()
        out_path.unlink()

    @staticmethod
    def test_bim_bgi_write():
        """Test writing bim.bgi"""
        bim_path = Path(Path(__file__).parent, "Data", "EUR.ldpred_21.bim")
        write_path = Path(Path(__file__).parent, "Data", "Write")
        PlinkObject(bim_path, bgi_write_path=write_path).create_bim_bgi()

        out_path = Path(Path(__file__).parent, "Data", "Write", "EUR.ldpred_21.bim.bgi")
        assert out_path.exists()
        out_path.unlink()

    def test_stats(self):
        """
        Check that we successfully parse the number of individuals and snps, check the length of arrays are equal to
        these stats values
        """
        bgen = BgenObject(Path(Path(__file__).parent, "Data", "EUR.ldpred_21.bgen"))

        self.assertEqual(bgen.iid_count, 483)
        self.assertEqual(bgen.sid_count, 7909)

        self.assertEqual(len(bgen.sid_array()), bgen.sid_count)
        self.assertEqual(len(bgen.iid_array()), bgen.iid_count)

    def test_parsers(self):
        """Test that loading the full array of snps actually returns all the snps"""
        bgen = BgenObject(Path(Path(__file__).parent, "Data", "EUR.ldpred_21.bgen"))

        self.assertEqual(bgen.sid_count, 7909)

        self.assertEqual(len(bgen.variant_array()), bgen.sid_count)
        self.assertEqual(len(bgen.info_array()), bgen.sid_count)
        self.assertEqual(len(bgen.dosage_array()), bgen.sid_count)

    def test_extractors(self):
        """Test all extractors work based on our three if elif else statments of len(snps) > 1, ==1, 0"""
        bgen = BgenObject(Path(Path(__file__).parent, "Data", "EUR.ldpred_21.bgen"))

        self.assertEqual(bgen.sid_count, 7909)

        snps = ['rs55776382', 'rs2801301']
        self.assertEqual(len(bgen.info_from_sid(snps)), 2)
        self.assertEqual(len(bgen.variant_from_sid(snps)), 2)
        self.assertEqual(len(bgen.dosage_from_sid(snps)), 2)

        # We originally had an issue where a single snp would cause a failure due to IN not being a valid sql call on a
        # tuple of length 1 hence validation like this
        snps = ['rs55776382']
        self.assertEqual(len(bgen.info_from_sid(snps)), 1)
        self.assertEqual(len(bgen.variant_from_sid(snps)), 1)
        self.assertEqual(len(bgen.dosage_from_sid(snps)), 1)

        snps = []
        self.assertEqual(len(bgen.info_from_sid(snps)), 0)
        self.assertEqual(len(bgen.variant_from_sid(snps)), 0)
        self.assertEqual(len(bgen.dosage_from_sid(snps)), 0)


if __name__ == '__main__':
    unittest.main()
