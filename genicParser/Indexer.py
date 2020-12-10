from pathlib import Path
import genicParser.errors as ec


class Bgi:
    def __init__(self, file_path):
        if Path(file_path).suffix == ".bgen":
            self.bgen = Path(file_path)
            self.bed_file = self.bim_file = self.fam_file = None
        else:
            self.bed_file, self.bim_file, self.fam_file, self.bgen = self.validate_paths(file_path)

    def a(self):
        assert self.bed_file

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
