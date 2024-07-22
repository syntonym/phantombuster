import zlib
import re
import numpy as np


class EndOfFile(Exception):
    pass

def read_header(f):
    header = f.read(12)
    if header == b"":
        raise EndOfFile()
    xlen = parse_xlen(header)

    subfields = f.read(xlen)
    if subfields == b"":
        raise EndOfFile()

    bsize = parse_bsize(subfields, xlen)
    return xlen, bsize

def read_block(f, bsize, xlen):
    """Assumes you called read_header before"""
    rest = f.read(bsize - 12 - xlen)
    return zlib.decompress(rest, -15)

def parse_xlen(header):
    return int.from_bytes(header[-2:], "little")

def parse_bsize(subfields, xlen):
    offset = 0
    block_size = None
    while offset <= xlen:
        if subfields[offset:offset+2] == b"BC":
            block_size = int.from_bytes(subfields[offset+4:offset+6], "little") +1
            break
        else:
            offset += int.from_bytes(subfields[offset+2:offset+4], "little") + 4
    if block_size is None:
        raise Exception("Could not find BSIZE")

    return block_size

def parse_read_start(data):
    if data == b"":
        raise EndOfFile()

    found = None
    i = 0
    for m in re.finditer(8*b"\xff", data):
        if found:
            offset = m.start() - found 
            if offset == 20 and (found-4 >= 0):
                return found-4 
            else:
                found = m.start()
        else:
            found = m.start()
        i += 1
        if i == 10:
            raise Exception("Could not find the beginning of the read, either read contains FF bytes or assumptions are violated")
    raise Exception("Could not find the beginning of the read, assumptions are violated")

def index(f, every):
    coffset = 0
    index = []
    global_i = 0
    i = 0
    first_block = True

    while True:
        i += 1
        global_i += 1
        f.seek(coffset)
        try:
            xlen, bsize = read_header(f)
        except EndOfFile:
            break

        #if True:
        if i % every == 0 or global_i <= 2:
            block_data = read_block(f, bsize, xlen)
            if first_block:
                uoffset = 0
            else:
                try:
                    uoffset = parse_read_start(block_data)
                except EndOfFile:
                    break
            idx = (coffset << 16) | uoffset
            if idx >= 0:
                index.append(idx)

        coffset += bsize
        first_block = False

    index = index[1:] + [np.inf]

    return index


