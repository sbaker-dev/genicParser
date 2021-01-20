
def path_invalid(parent_path, operation):
    return f"INVALID Path for {operation}\n" \
           f"{operation} attempt to navigate to a directory or file at {parent_path} but it does not exist"


def magic_violation(file_name):
    return f"INVALID BGEN FILE for file at path: {file_name}\n" \
           f"Bgen files have a magic number, which in more modern files is 4 bytes of b'bgen' but can" \
           f" also take four zero bytes in older files.\n"


def offset_violation(file_name, offset, header_block_length):
    return f"INVALID HEADER BLOCK for file at path: {file_name}\n" \
           f"Bgen files have a header block, that must not be larger than the offset. Yet found\n" \
           f"Offset: {offset} header_block_length {header_block_length}"


def sample_block_violation(header, offset, block_size):
    return f"INVALID BLOCK SIZE\n" \
           f"The header block + the offset should equal the length of the sample block yet found\n" \
           f"header {header}, offset {offset}, block_size {block_size}"


def sample_size_violation(sample_size, found_size):
    return f"INVALID SAMPLE SIZE\n" \
           f"The sample size for this file is set to {sample_size} yet the file found {found_size}"


def compression_violation(file_name, compression_flag):
    return f"INVALID COMPRESSION FLAG for file at path: {file_name}\n" \
           f"Bgen files have a flag, where the first two bits represent the compression of the data with 0 being " \
           f"uncompressed, 1 being compressed via zlib, and 2 being compressed via z-standard. Yet found\n" \
           f"Compression flag: {compression_flag}"


def layout_violation(file_name, layout_flag):
    return f"INVALID LAYOUT FLAG for file at path: {file_name}\n" \
           f"Bgen files have a flag, where bits 2-5 represents a Layout Flags should only take a value of 1 or 2 and " \
           f"relates to how and what genetic information is stored. Yet Found\n" \
           f"Layout flag: {layout_flag}"


def sample_identifier_violation(file_name, sample_identifier):
    return f"INVALID SAMPLE IDENTIFIER FLAG for file at path: {file_name}\n" \
           f"Bgen files have a flag, where bit 31 represents if sample identifiers as within the file at 1 or not" \
           f" at 0. Yet found\n" \
           f"Sample Identifier flag: {sample_identifier}"


def bgi_path_violation(bgi_path):
    return f"INVALID BGI TYPE for bgi_path of type: {bgi_path}\n" \
           f"Bgi_path defaults to False where you don't have an index in an external file. If you do have .bgi file" \
           f" in the same directory, name the same thing but with .bgi on the end, then it should be True. If the " \
           f"file is in another directory you may also pass the path as a str in to bgi_path. Yet Found\n" \
           f"BGI path type: {type(bgi_path)}"


def index_violation(operation):
    return f"NO INDEX HAS BEEN SET\n" \
           f"Operation {operation} requires you to have a .bgi file which acts as a sql database"


def bgen_no_variant_found(variant_name):
    return f"NO VARIANT NAMED {variant_name} FOUND"


def invalid_slice(slice_item):
    return f"INVALID SLICE\n" \
           f"Slice takes a slice or a tuple of two slices yet was passed type - {type(slice_item)} for {slice_item}"


def slice_error(slice_type, slice_length):
    return f"TUPLE OF SLICES IS NOT EQUAL TO 2\n" \
           f"You must pass two slices for iid and sid, for example [:7, :7] for the first 7 iid and snps\n" \
           f"However found a {slice_type} of length {slice_length}"


def wrong_slice_type(slice_type):
    return f"SLICE IS NETHER SLICE OR LIST\n" \
           f"You must pass the indexer of iid or sid as a slice, for example :3, or as a list of indexes, eg[0, 1, 2]" \
           f"\nYet Found {slice_type}"


def slice_list_type():
    return f"SLICE LIST IS NOT A LIST OF INTS\n" \
           f"Slicing takes a slice, or a list of ints that are the indexes.\nDid you forget to convert iids/sids into " \
           f"indexes?"
