def bits_to_int(bits):
    """Converts bits to int."""
    result = 0
    for bit in bits:
        result = (result << 1) | bit

    return result


def no_decompress(data):
    """Don't decompress"""
    return data
