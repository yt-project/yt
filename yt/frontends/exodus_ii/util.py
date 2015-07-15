import string
from itertools import takewhile
from collections import OrderedDict
from netCDF4 import Dataset
import numpy as np
import re

class ExodusIIData:
    def __init__(self, filename):
        self.filename = filename
        self.dataset = OrderedDict()
        self.coordinates = np.array([])
        self.connectivity = []
        self.info_records = {}
        self.data = []

    def read(self):
        print "Reading netCDF4 Dataset"
        self.dataset = Dataset(self.filename).variables

        print "Loading coordinates"
        self.coordinates = np.array([self.dataset["coord%s" % ax][:]
                                     for ax in 'xyz']).transpose().copy()

        print "Loading connectivity"
        for i in range(self.get_num_elem()):
            self.add_connect("connect%s" % (i+1))
            ci = self.connectivity[-1]
            vals = {}

            for j, var_name in enumerate(self.get_var_names()):
                vals['gas', var_name] = self.dataset["vals_elem_var%seb%s" % (j+1, i+1)][:].astype("f8")[-1,:]

            for j, nod_name in enumerate(self.get_nod_names()):
                # We want just for this set of nodes all the node variables
                # Use (ci - 1) to get these values
                vals['gas', nod_name] = self.dataset["vals_nod_var%s" % (j+1)][:].astype("f8")[-1, ci - 1, ...]

            self.data.append(vals)

        if 'info_records' in self.dataset.keys():
            print "Loading info records"
            self.info_records = load_info_records(self.dataset['info_records'])

    def get_num_elem(self):
        return self.dataset["eb_status"][:].shape[0]

    def add_connect(self, connect_key):
        self.connectivity.append(self.dataset[connect_key][:].astype("i8"))

    def get_var_names(self):
        return [sanitize_string(v.tostring()) for v in
                self.dataset["name_elem_var"][:]]

    def get_nod_names(self):
        return [sanitize_string(v.tostring()) for v in
                self.dataset["name_nod_var"][:]]



def sanitize_string(s):
    return "".join(_ for _ in takewhile(lambda a: a in string.printable, s))

def load_info_records(info_records):
    info_records_parsed = ["".join(line_chars) for line_chars in info_records]
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

    version_info = OrderedDict()

    for line in grouped['Version Info']:
        split_line = line.split(":")
        key = split_line[0]
        val = ":".join(split_line[1:]).lstrip().rstrip()
        if key != '':
            version_info[key] = val

    grouped['Version Info'] = version_info
    
    return grouped
    
def get_top_levels(info_records):
    top_levels = []
    for idx, line in enumerate(info_records):
        pattern = re.compile("###[a-zA-Z\s]+")
        if pattern.match(line):
            clean_line = line.strip('###').lstrip().rstrip()
            top_levels.append([idx, clean_line])
    
    return top_levels


