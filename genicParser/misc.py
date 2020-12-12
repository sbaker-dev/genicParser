from . import errors_codes as ec
from math import ceil
import numpy as np
import struct
import os


def set_bgi(bgi_present, base_file_path):
    """
    Connect to the index file either via a bgi file in the same directory or in another directory.
    """

    if not bgi_present:
        return False
    elif isinstance(bgi_present, str):
        if not os.path.isfile(f"{bgi_present}"):
            raise IOError(f"{bgi_present} was not found")
        else:
            return True
    elif bgi_present:
        if not os.path.isfile(f"{base_file_path}.bgi"):
            raise IOError(f"{base_file_path}.bgi was not found")
        else:
            return True
    else:
        raise TypeError(ec.bgi_path_violation(bgi_present))


def bits_to_int(bits):
    """Converts bits to int."""
    result = 0
    for bit in bits:
        result = (result << 1) | bit

    return result


def byte_to_int(byte):
    """Converts a byte to a int for python 3."""
    return byte


def pack_bits(data, b):
    """Unpacks BGEN probabilities (as bits)."""
    # Getting the data from the bytes
    data = np.fromiter(
        ((byte_to_int(byte) >> i) & 1 for byte in data for i in range(8)),
        dtype=bool,
    )
    data.shape = (data.shape[0] // b, b)

    # Finding the closest full bytes (if required)
    # TODO: Improve this so that it is more efficient
    full_bytes = data[:, ::-1]
    if data.shape[1] % 8 != 0:
        nb_bits = int(ceil(b / 8)) * 8
        full_bytes = np.zeros((data.shape[0], nb_bits), dtype=bool)
        full_bytes[:, -b:] += data[:, ::-1]

    # Packing the bits
    packed = np.packbits(full_bytes, axis=1)

    # Left-shifting for final value
    final = packed[:, 0]
    for i in range(1, packed.shape[1]):
        final = np.left_shift(final, 8, dtype=np.uint) | packed[:, i]

    return final


def struct_unpack(struct_format, data, list_return=False):
    if list_return:
        return struct.unpack(struct_format, data)
    else:
        return struct.unpack(struct_format, data)[0]


def no_decompress(data):
    """Don't decompress"""
    return data
