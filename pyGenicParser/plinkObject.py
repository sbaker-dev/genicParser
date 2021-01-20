from .variantObjects import BimVariant, FamId, Variant
from . import errors_codes as ec
from . import misc as mc

from pathlib import Path
import numpy as np
import sqlite3


class PlinkObject:
    def __init__(self, genetic_path, bgi_present=False):
        self.bed_file_path, self.bim_file_path, self.fam_file_path = self.validate_paths(genetic_path)
        self.bim_file = open(self.bim_file_path, "r")
        self.fam_file = open(self.fam_file_path, "r")

        # Set the bgi file if present, and store this for indexing if required.
        self.bgi_present = bgi_present
        self.bgi_file = mc.set_bgi(self.bgi_present, self.bim_file_path)
        if self.bgi_file:
            self.bim_connection, self.bim_index = self._connect_to_bgi_index()
        else:
            self.bim_connection, self.bim_index = None, None

    def close_all(self):
        """Close all open files"""
        self.bim_file.close()
        self.fam_file.close()

    def info_array(self, as_variant=False):
        """Return an array of all the variants in the bgen file"""
        assert self.bim_index, ec.index_violation("info_array")
        if as_variant:
            self.bim_index.execute("SELECT chromosome, position, rsid, allele1, allele2 FROM Variant")
            return np.array([Variant(chromosome, position, snp_id, a1, a2) for chromosome, position, snp_id, a1, a2
                             in self.bim_index.fetchall()])
        else:
            self.bim_index.execute("SELECT chromosome, rsid, morgan_pos, position, allele1, allele2 FROM Variant")
            return np.array([BimVariant(chromosome, variant_id, morgan_pos, bp_position, a1, a2)
                             for chromosome, variant_id, morgan_pos, bp_position, a1, a2 in self.bim_index.fetchall()])

    def info_from_sid(self, snp_names, as_variant=False):
        """Construct an array of variant identifiers for all the snps provided to snp_names"""
        assert self.bim_index, ec.index_violation("variant_info_from_sid")
        if as_variant:
            self.bim_index.execute("SELECT chromosome, position, rsid, allele1, allele2 FROM Variant"
                                   " WHERE rsid IN {}".format(tuple(snp_names)))
            return np.array([Variant(chromosome, position, snp_id, a1, a2) for chromosome, position, snp_id, a1, a2
                             in self.bim_index.fetchall()])
        else:
            self.bim_index.execute("SELECT chromosome, rsid, morgan_pos, position, allele1, allele2 FROM Variant"
                                   " WHERE rsid IN {}".format(tuple(snp_names)))
            return np.array([BimVariant(chromosome, variant_id, morgan_pos, bp_position, a1, a2)
                             for chromosome, variant_id, morgan_pos, bp_position, a1, a2 in self.bim_index.fetchall()])

    def create_bim_bgi(self, bgi_write_path=None):
        """
        This will create a 'mock' .bgi akin to bgenix but with a few differences. Firstly, given information of plink is
        stored in different files this new .bgi acts as the old .bim. It contains all the information bim does, but with
        the variant starting position within the bed file so that it can quickly be parsed out.

        This also contains some misc data such as the count of iid and sid so that it can be quickly accessed.
        """
        # Check if the file already exists
        if not bgi_write_path:
            write_path = Path(str(self.bim_file_path.absolute()) + ".bgi")
        else:
            write_path = str(Path(bgi_write_path, self.bim_file_path.name).absolute()) + ".bgi"

        if Path(write_path).exists():
            print(f"Bgi Already exists for {self.bim_file_path.name}")
            self.close_all()
        else:
            # Establish the connection
            connection = sqlite3.connect(write_path)
            c = connection.cursor()

            # Load the bim data as a bgi index
            bim_dict = self.construct_bim_index(bgi_index=True)

            # Set the number of snps to the be the length of the dict, and get the length of iid from fam length
            sid_count = len(bim_dict)
            iid_count = len(self.get_family_identifiers())
            self.close_all()

            # Construct the bed array based on its byte formula
            # See https://www.cog-genomics.org/plink/1.9/formats#bed
            bed_array = [int(np.ceil(0.25 * iid_count) * bimIndex + 3) for bimIndex in np.arange(sid_count)]

            # Append this to the front of our dict of values
            for index, k in enumerate(bim_dict):
                bim_dict[k] = [bed_array[index]] + bim_dict[k]

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

    def _connect_to_bgi_index(self):
        """Connect to the index (which is an SQLITE database)."""
        bim_file = sqlite3.connect(str(self.bim_file_path.absolute()) + ".bgi")
        return bim_file, bim_file.cursor()

    def construct_bim_index(self, bgi_index=False):
        """
        Bim files need to be index via seek, so we can extract a given snp loci without having to store all of this of
        then in memory
        """
        indexer = {}
        cumulative_seek = 0
        for line in self.bim_file:
            chromosome, variant_id, morgan_pos, bp_position, a1, a2 = line.split()
            if bgi_index:
                indexer[variant_id] = [cumulative_seek, variant_id, chromosome, morgan_pos, bp_position, a1, a2]
            else:
                indexer[variant_id] = cumulative_seek
            cumulative_seek += len(line)

        return indexer

    def get_variant(self, seek, as_variant=False):
        """
        Extract a given variants loci based on the seek index from construct_index
        :param seek: The seek index
        :param as_variant: If you want it as a standardised across parameter variant, or a Bim Variant with morgan pos
        :return: The line
        """
        self.bim_file.seek(seek)
        chromosome, variant_id, morgan_pos, bp_position, a1, a2 = self.bim_file.readline().split()

        if as_variant:
            return BimVariant(chromosome, variant_id, morgan_pos, bp_position, a1, a2).to_variant()
        else:
            return BimVariant(chromosome, variant_id, morgan_pos, bp_position, a1, a2)

    def get_family_identifiers(self):
        """
        This will iterate through the fam file and extract the information
        """
        fam_data = []
        for line in self.fam_file:
            fid, iid, i_fid, i_mid, sex, phenotype = line.split()
            fam_data.append(FamId(fid, iid, i_fid, i_mid, sex, phenotype))

        return fam_data

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

        return bed, bim, fam
