import base64
import re
import string
import zlib
from collections import OrderedDict
from itertools import takewhile

import numpy as np

type_decider = {"Float32": "<f4", "Float64": "<f8"}


def decode_piece(xmlPiece):
    coord_type = type_decider[xmlPiece["Points"]["DataArray"]["@type"]]

    _, coords = decode_binary(
        xmlPiece["Points"]["DataArray"]["#text"].encode(), dtype=coord_type
    )
    _, conn = decode_binary(
        xmlPiece["Cells"]["DataArray"][0]["#text"].encode(), dtype="u4"
    )
    _, offsets = decode_binary(
        xmlPiece["Cells"]["DataArray"][1]["#text"].encode(), dtype="u4"
    )
    _, cell_types = decode_binary(
        xmlPiece["Cells"]["DataArray"][2]["#text"].encode(), dtype="u1"
    )

    coords = coords.reshape((coords.size // 3, 3))

    return coords, conn, offsets, cell_types


def decode_binary(blob, use_zlib=True, dtype="<f4"):
    split_location = blob.find(b"==") + 2
    first = base64.decodebytes(blob[:split_location])
    second = base64.decodebytes(blob[split_location:])
    if zlib:
        second = zlib.decompress(second)
    return np.frombuffer(first, dtype="<f4"), np.frombuffer(second, dtype=dtype)


def get_num_pseudo_dims(coords):
    D = coords.shape[1]
    return sum([np.all(coords[:, dim] == 0.0) for dim in range(D)])


def sanitize_string(s):
    _printable = {ord(_) for _ in string.printable}
    return "".join([chr(_) for _ in takewhile(lambda a: a in _printable, s)])


def load_info_records(info_records):
    info_records_parsed = [sanitize_string(line_chars) for line_chars in info_records]
    return group_by_sections(info_records_parsed)


def group_by_sections(info_records):
    # 1. Split by top groupings
    top_levels = get_top_levels(info_records)
    # 2. Determine if in section by index number
    grouped = OrderedDict()
    for tidx, top_level in enumerate(top_levels):
        grouped[top_level[1]] = []

        try:
            next_idx = top_levels[tidx + 1][0]
        except IndexError:
            next_idx = len(info_records) - 1

        for idx in range(top_level[0], next_idx):
            if idx == top_level[0]:
                continue

            grouped[top_level[1]].append(info_records[idx])

    if "Version Info" in grouped.keys():
        version_info = OrderedDict()
        for line in grouped["Version Info"]:
            split_line = line.split(":")
            key = split_line[0]
            val = ":".join(split_line[1:]).lstrip().rstrip()
            if key != "":
                version_info[key] = val
        grouped["Version Info"] = version_info

    return grouped


def get_top_levels(info_records):
    top_levels = []
    for idx, line in enumerate(info_records):
        pattern = re.compile(r"###[a-zA-Z\s]+")
        if pattern.match(line):
            clean_line = re.sub(r"[^\w\s]", "", line).lstrip().rstrip()
            top_levels.append([idx, clean_line])

    return top_levels
