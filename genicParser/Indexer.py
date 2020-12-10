from pathlib import Path
import genicParser.errors as ec
import sqlite3
import numpy as np


class Bgi:
    def __init__(self, file_path):
        if Path(file_path).suffix == ".bgen":
            self.bgen = Path(file_path)
            self.bed_file = self.bim_file = self.fam_file = None
        else:
            self.bed_file, self.bim_file, self.fam_file, self.bgen = self.validate_paths(file_path)

    def a(self):
        assert self.bed_file
        connection = sqlite3.connect(r"C:\Users\Samuel\Documents\Genetic_Examples\PolyTutOut\BgenLoader\test.bgi")
        c = connection.cursor()

        bim_dict = {}
        with open(self.bim_file, "r") as f:
            cumulative_seek = 0
            for line in f:
                chromosome, variant_id, morgan_pos, bp_position, a1, a2 = line.split()
                bim_dict[variant_id] = [cumulative_seek, variant_id, chromosome, morgan_pos, bp_position, a1, a2]

                cumulative_seek += len(line)

        sid_count = len(bim_dict)
        iid_count = len([la for la in open(self.fam_file, "r")])

        arrayaaa = [int(np.ceil(0.25 * iid_count) * bimIndex + 3) for bimIndex in np.arange(sid_count)]

        for index, k in enumerate(bim_dict):
            bim_dict[k] = [arrayaaa[index]] + bim_dict[k]

        print(bim_dict["rs9425291"])

        c.execute('''CREATE TABLE Variant (
            bed_start_position INTEGER,
            bim_start_position INTEGER,
            rsid TEXT,
            chromosome INTEGER,
            morgan_pos REAL,
            position INTEGER,
            allele1 TEXT,
            allele2 TEXT
        )''')

        for value in bim_dict.values():

            c.execute(f'INSERT INTO Variant VALUES {tuple(value)}')

        c.execute('''CREATE TABLE Misc (
            sid_count INTEGER,
            iid_count INTEGER
        )''')

        c.execute(f'INSERT INTO Misc VALUES {tuple([sid_count, iid_count])}')

        connection.commit()
        connection.close()


    def b(self):
        connection = sqlite3.connect(r"C:\Users\Samuel\Documents\Genetic_Examples\PolyTutOut\BgenLoader\test.bgi")
        c = connection.cursor()

        c.execute("SELECT bed_start_position, rsid FROM Variant")
        available_table = (c.fetchall())
        print(available_table)

        c.execute("SELECT sid_count FROM Misc")

        available_table = (c.fetchall())
        print(available_table)


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
