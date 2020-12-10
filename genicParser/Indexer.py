from .plinkObject import PlinkObject
from .bgenObject import BgenObject

from pathlib import Path
import numpy as np
import sqlite3


class Bgi:
    def __init__(self, file_path):
        self.file_path = Path(file_path)

    def create_bim_bgi(self):
        """
        This will create a 'mock' .bgi akin to bgenix but with a few differences. Firstly, given information of plink is
        stored in different files this new .bgi acts as the old .bim. It contains all the information bim does, but with
        the variant starting position within the bed file so that it can quickly be parsed out.

        This also contains some misc data such as the count of iid and sid so that it can be quickly accessed.
        """
        assert self.file_path.suffix == (".bed" or ".bim" or ".fam")

        # Construct a plink object holder
        plink = PlinkObject(self.file_path)

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
        assert self.file_path.suffix == ".bgen"

        # Can only work on bgen 1.2, so validate this is true
        bgen_object = BgenObject(self.file_path)
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
