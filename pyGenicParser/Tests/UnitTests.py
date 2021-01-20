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
    def _loader():
        """Call the bgen file"""
        return BgenObject(Path(Path(__file__).parent, "Data", "EUR.ldpred_21.bgen"))

    def test_bgen_bgi_write(self):
        """Test writing bgen.bgi"""
        write_path = Path(Path(__file__).parent, "Data", "Write")
        self._loader().create_bgi(write_path)

        out_path = Path(Path(__file__).parent, "Data", "Write", "EUR.ldpred_21.bgen.bgi")
        assert out_path.exists()
        out_path.unlink()

    @staticmethod
    def test_bim_bgi_write():
        """Test writing bim.bgi"""
        bim_path = Path(Path(__file__).parent, "Data", "EUR.ldpred_21.bim")
        write_path = Path(Path(__file__).parent, "Data", "Write")
        PlinkObject(bim_path, write_path).create_bim_bgi(write_path)

        out_path = Path(Path(__file__).parent, "Data", "Write", "EUR.ldpred_21.bim.bgi")
        assert out_path.exists()
        out_path.unlink()

    def test_stats(self):
        """
        Check that we successfully parse the number of individuals and snps, check the length of arrays are equal to
        these stats values
        """
        bgen = self._loader()

        self.assertEqual(bgen.iid_count, 483)
        self.assertEqual(bgen.sid_count, 7909)

        self.assertEqual(len(bgen.sid_array()), bgen.sid_count)
        self.assertEqual(len(bgen.iid_array()), bgen.iid_count)

    def test_parsers(self):
        """Test that loading the full array of snps actually returns all the snps"""
        bgen = self._loader()

        self.assertEqual(bgen.sid_count, 7909)

        self.assertEqual(len(bgen.variant_array()), bgen.sid_count)
        self.assertEqual(len(bgen.info_array()), bgen.sid_count)
        self.assertEqual(len(bgen.dosage_array()), bgen.sid_count)

    def test_extractors(self):
        """Test all extractors work based on our three if elif else statements of len(snps) > 1, ==1, 0"""
        bgen = self._loader()
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
        print("Should return 3 no names passed")
        self.assertEqual(len(bgen.info_from_sid(snps)), 0)
        self.assertEqual(len(bgen.variant_from_sid(snps)), 0)
        self.assertEqual(len(bgen.dosage_from_sid(snps)), 0)

    def test_indexing(self):
        """Tests if indexing on snp or iid actually produces an indexed request"""
        bgen = self._loader()
        sid_array = ['rs55776382', 'rs2801301', 'rs3869758', 'rs79913394', 'rs35829851', 'rs118189563']

        self.assertEqual(bgen[:, [1, 2, 3]].sid_count, 3)
        self.assertEqual(bgen[[1, 2, 3], :].iid_count, 3)
        self.assertEqual(bgen[:, :3].sid_count, 3)
        self.assertEqual(bgen[:3, :].iid_count, 3)
        self.assertEqual(bgen[:, bgen.sid_to_index(sid_array)].sid_count, 6)

        # Check we still get the right length even with return false (returns -1 if we fail to find something)
        bgen = bgen[:, bgen.sid_to_index(sid_array[:3])]
        self.assertEqual(bgen[:, bgen.sid_to_index(sid_array, True)].sid_count, 3)

        # Now do the same with iid
        iid_list = [[0, 0], [1, 1], [2, 2]]
        bgen = bgen[bgen.iid_to_index(iid_list), :]
        self.assertEqual(bgen.iid_count, 3)

    def test_indexed_data(self):

        bgen = self._loader()
        sid_array = ['rs55776382', 'rs2801301', 'rs3869758', 'rs79913394', 'rs35829851', 'rs118189563']

        bgen = bgen[:, bgen.sid_to_index(sid_array[:3])]

        self.assertEqual(len(bgen.sid_array()), 3)
        self.assertEqual(len(bgen.info_array()), 3)
        self.assertEqual(len(bgen.dosage_array()), 3)
        self.assertEqual(len(bgen.variant_array()), 3)

        # Check that snp names > than length of bgen.sid_count don't cause a crash
        self.assertEqual(len(bgen.info_from_sid(sid_array)), 3)
        self.assertEqual(len(bgen.dosage_from_sid(sid_array)), 3)
        self.assertEqual(len(bgen.variant_from_sid(sid_array)), 3)


if __name__ == '__main__':
    unittest.main()
