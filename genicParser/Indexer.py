from .plinkObject import PlinkObject
from .bgenObject import BgenObject
from . import errors_codes as ec

from pathlib import Path
import numpy as np
import sqlite3


class Bgi:
    def __init__(self, file_path):
        if Path(file_path).suffix == ".bgen":
            self.bgen = Path(file_path)
            self.bed_file = self.bim_file = self.fam_file = None
        else:
            self.bed_file, self.bim_file, self.fam_file, self.bgen = self.validate_paths(file_path)

    def create_bim_bgi(self):
        """
        This will create a 'mock' .bgi akin to bgenix but with a few differences. Firstly, given information of plink is
        stored in different files this new .bgi acts as the old .bim. It contains all the information bim does, but with
        the variant starting position within the bed file so that it can quickly be parsed out.

        This also contains some misc data such as the count of iid and sid so that it can be quickly accessed.
        """
        assert self.bed_file

        # Construct a plink object holder
        plink = PlinkObject(self.bim_file)

        # Load the bim data as a bgi index
        bim_dict = plink.construct_bim_index(bgi_index=True)

        # Set the number of snps to the be the length of the dict, and get the length of iid from fam length
        sid_count = len(bim_dict)
        iid_count = len(plink.get_family_identifiers())
        plink.close_all()

        # Construct the bed array based on its btye formula
        # See https://www.cog-genomics.org/plink/1.9/formats#bed
        bed_array = [int(np.ceil(0.25 * iid_count) * bimIndex + 3) for bimIndex in np.arange(sid_count)]

        # Append this to the front of our dict of values
        for index, k in enumerate(bim_dict):
            bim_dict[k] = [bed_array[index]] + bim_dict[k]

        # Establish the connection
        connection = sqlite3.connect(r"C:\Users\Samuel\Documents\Genetic_Examples\PolyTutOut\BgenLoader\test.bgi")
        c = connection.cursor()

        # Create our core table that mimics bgi from bgenix but with bed and bim
        c.execute('''
            CREATE TABLE Variant (
            bed_start_position INTEGER,
            bim_start_position INTEGER,
            rsid TEXT,
            chromosome INTEGER,
            morgan_pos REAL,
            position INTEGER,
            allele1 TEXT,
            allele2 TEXT
        )''')

        # Append our values into this table
        for value in bim_dict.values():
            c.execute(f'INSERT INTO Variant VALUES {tuple(value)}')

        # Create a misc table of sid_count and iid_count
        c.execute('''
            CREATE TABLE Misc (
            sid_count INTEGER,
            iid_count INTEGER
        )''')
        c.execute(f'INSERT INTO Misc VALUES {tuple([sid_count, iid_count])}')

        # Commit the file
        connection.commit()
        connection.close()

    def create_bgen_bgi(self):
        """
        THis attempts to directly mimic bgi but without the meta data table so that these files can be create without
        bgenix, but also allows for it to be cross compatible with files that have been made via bgenix
        """
        assert self.bgen

        # Can only work on bgen 1.2, so validate this is true
        bgen_object = BgenObject(self.bgen)
        assert bgen_object.layout == 2

        # Create our lines for the table
        bgi_lines = bgen_object.create_bgi()

        # Establish the connection
        connection = sqlite3.connect(r"C:\Users\Samuel\Documents\Genetic_Examples\PolyTutOut\BgenLoader\test2.bgi")
        c = connection.cursor()

        # Create our core table that mimics Variant bgi from bgenix
        c.execute('''
            CREATE TABLE Variant (
            file_start_position INTEGER,
            size_in_bytes INTEGER,
            chromosome INTEGER,
            position INTEGER,
            rsid TEXT,
            allele1 TEXT,
            allele2 TEXT
              )''')

        # Write values to table
        for value in bgi_lines:
            c.execute(f'INSERT INTO Variant VALUES {tuple(value)}')

        # Commit the file
        connection.commit()
        connection.close()

    @staticmethod
    def validate_paths(genetic_path):
        """
        Users may submit a path to a specific file within plink, such as a .bed/.bim/.fam or they just provide the root
        name. This method validates and returns the paths.

        :return: The path to the bed, bim, and fam file in that order
        """
        # Construct path as an object
        ld_path = Path(genetic_path)

        # Check file home directory can be reached
        assert ld_path.parent.exists(), ec.path_invalid(ld_path.parent, "_set_ld_ref")

        # If a file has a plink suffix take the stem of the name otherwise just take the name
        if (ld_path.suffix == ".bed") or (ld_path.suffix == ".bim") or (ld_path == ".fam"):
            bed = Path(f"{str(ld_path.parent)}/{ld_path.stem}.bed")
            bim = Path(f"{str(ld_path.parent)}/{ld_path.stem}.bim")
            fam = Path(f"{str(ld_path.parent)}/{ld_path.stem}.fam")
        else:
            bed = Path(f"{str(ld_path.parent)}/{ld_path.name}.bed")
            bim = Path(f"{str(ld_path.parent)}/{ld_path.name}.bim")
            fam = Path(f"{str(ld_path.parent)}/{ld_path.name}.fam")

        # Check the files exists then return with mode of plink, no bgen object and a bed, bim and fam file.
        assert bed.exists(), ec.path_invalid(bed, "_set_ld_ref")
        assert bim.exists(), ec.path_invalid(bim, "_set_ld_ref")
        assert fam.exists(), ec.path_invalid(fam, "_set_ld_ref")

        return bed, bim, fam, None
