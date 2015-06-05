import string
from itertools import takewhile
from collections import OrderedDict
import numpy as np
from netCDF4 import Dataset

class ExodusData:
    def __init__(self, filename):
        self.filename = filename
        self.dataset = OrderedDict()
        self.coords = np.array([])
        self.connects = []
        self.data = []

    def load_dataset(self):
        self.dataset = Dataset(self.filename).variables

        self.coords = np.array([self.dataset["coord%s" % ax][:]
                                for ax in 'xyz']).transpose().copy()

        for i in range(self.get_num_elem()):
            self.add_connect("connect%s" % (i+1))
            ci = self.connects[-1]
            vals = {}

            for j, var_name in enumerate(self.get_var_names()):
                vals['gas', var_name] = self.dataset["vals_elem_var%seb%s" % (j+1, i+1)][:].astype("f8")[-1,:]

            for j, nod_name in enumerate(self.get_nod_names()):
                # We want just for this set of nodes all the node variables
                # Use (ci - 1) to get these values
                vals['gas', nod_name] = self.dataset["vals_nod_var%s" % (j+1)][:].astype("f8")[-1, ci - 1, ...]

            self.data.append(vals)

    def get_num_elem(self):
        return self.dataset["eb_status"][:].shape[0]

    def add_connect(self, connect_key):
        self.connects.append(self.dataset[connect_key][:].astype("i8"))

    def get_var_names(self):
        return [self.sanitize_string(v.tostring()) for v in
                self.dataset["name_elem_var"][:]]

    def get_nod_names(self):
        return [self.sanitize_string(v.tostring()) for v in
                self.dataset["name_nod_var"][:]]

    def sanitize_string(self, s):
        return "".join(_ for _ in takewhile(lambda a: a in string.printable, s))

