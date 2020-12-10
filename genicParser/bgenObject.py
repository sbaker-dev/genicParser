from . import misc as mc
from . import errors_codes as ec
from pathlib import Path
import numpy as np
import struct
import zlib
import zstd


class BgenObject:
    def __init__(self, file_path, bgi_present=True, probability=None, iter_array_size=1000):
        """

        :param file_path:

        :param bgi_present: Takes a value of True if the .bgi is in the same directory and named file_path.bgi
            otherwise can ec passed as a path if it is in a different directory.
        :type bgi_present: bool | str

        :param probability:
        """

        self.file_path = Path(file_path)
        self._bgen_binary = open(file_path, "rb")

        self.offset, self.headers, self.sid_count, self.iid_count, self.compression, self.layout, \
            self.sample_identifiers, self._variant_start = self.parse_header()

        self.probability = probability
        self.iter_array_size = iter_array_size

        self.bgi_present = bgi_present

    def parse_header(self):
        """
        Extract information from the header of the bgen file.

        Spec at https://www.well.ox.ac.uk/~gav/bgen_format/spec/latest.html

        :return: offset, headers, variant_number, sample_number, compression, layout, and sample_identifiers
        """

        # Check the header block is not larger than offset
        offset = self.unpack("<I", 4)
        headers = self.unpack("<I", 4)
        assert headers <= offset, ec.offset_violation(self._bgen_binary.name, offset, headers)
        variant_start = offset + 4

        # Extract the number of variants and samples
        variant_number = self.unpack("<I", 4)
        sample_number = self.unpack("<I", 4)

        # Check the file is valid
        magic = self.unpack("4s", 4)
        assert (magic == b'bgen') or (struct.unpack("<I", magic)[0] == 0), ec.magic_violation(self._bgen_binary.name)

        # Skip the free data area
        self._bgen_binary.read(headers - 20)

        # Extract the flag, then set compression layout and sample identifiers from it
        compression, layout, sample_identifiers = self._header_flag()
        return offset, headers, variant_number, sample_number, compression, layout, sample_identifiers, variant_start

    def _header_flag(self):
        """
        The flag represents a 4 byte unsigned int, where the bits relates to the compressedSNPBlock at bit 0-1, Layout
        at 2-5, and sampleIdentifiers at 31

        Spec at https://www.well.ox.ac.uk/~gav/bgen_format/spec/latest.html

        :return: Compression, layout, sampleIdentifiers
        """
        # Reading the flag
        flag = np.frombuffer(self._bgen_binary.read(4), dtype=np.uint8)
        flag = np.unpackbits(flag.reshape(1, flag.shape[0]), bitorder="little")

        # [N1] Bytes are stored right to left hence the reverse, see shorturl.at/cOU78
        # Check the compression of the data
        compression_flag = mc.bits_to_int(flag[0: 2][::-1])
        assert 0 <= compression_flag < 3, ec.compression_violation(self._bgen_binary.name, compression_flag)
        if compression_flag == 0:
            compression = mc.no_decompress
        elif compression_flag == 1:
            compression = zlib.decompress
        else:
            compression = zstd.decompress

        # Check the layout is either 1 or 2, see [N1]
        layout = mc.bits_to_int(flag[2:6][::-1])
        assert 1 <= layout < 3, ec.layout_violation(self._bgen_binary.name, layout)

        # Check if the sample identifiers are in the file or not, then return
        assert flag[31] == 0 or flag[31] == 1, ec.sample_identifier_violation(self._bgen_binary.name, flag[31])
        if flag[31] == 0:
            return compression, layout, False
        else:
            return compression, layout, True

    def unpack(self, struct_format, size, list_return=False):
        """
        Use a given struct formatting to unpack a byte code

        Struct formatting
        ------------------
        https://docs.python.org/3/library/struct.html

        :param struct_format: The string representation of the format to use in struct format. See struct formatting for
            a list of example codes.
        :type struct_format: str

        :param size: The byte size
        :type size: int

        :key list_return: If we expect multiple values then we return a tuple of what was unpacked, however if there is
            only one element then we often just index the first element to return it directly. Defaults to false.
        :type list_return: bool

        :return: Whatever was unpacked
        :rtype: Any
        """
        if list_return:
            return struct.unpack(struct_format, self._bgen_binary.read(size))
        else:
            return struct.unpack(struct_format, self._bgen_binary.read(size))[0]
