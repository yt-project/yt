import re
import string
from collections import OrderedDict
from itertools import takewhile

import numpy as np


def get_num_pseudo_dims(coords):
    D = coords.shape[1]
    return sum(np.all(coords[:, dim] == 0.0) for dim in range(D))


def sanitize_string(s):
    _printable = {ord(_) for _ in string.printable}
    return "".join(chr(_) for _ in takewhile(lambda a: a in _printable, s))


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
