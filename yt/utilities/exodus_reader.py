import string
from itertools import takewhile
import numpy as np
from netCDF4 import Dataset

class ExodusData:
    def __init__(self, filename):
        self.filename = filename
        self.dataset = []
        self.coords = []
        self.connects = []
        self.data = []

    def load_dataset(self):
        self.dataset = Dataset(self.filename).variables
        nelem = self.dataset["eb_status"][:].shape[0]
        coord = np.array([self.dataset["coord%s" % ax][:]
                         for ax in 'xyz']).transpose().copy()

        for i in range(nelem):
            self.connects.append(self.dataset["connect%s" % (i+1)][:].astype("i8"))
            ci = self.connects[-1]
            self.coords.append(coord)  # Same for all
            vals = {}

            for j, v in enumerate(self.get_var_names()):
                vals['gas', v] = self.dataset["vals_elem_var%seb%s" % (j+1, i+1)][:].astype("f8")[-1,:]

            for j, v in enumerate(self.get_nod_names()):
                # We want just for this set of nodes all the node variables
                # Use (ci - 1) to get these values
                vals['gas', v] = self.dataset["vals_nod_var%s" % (j+1)][:].astype("f8")[-1, ci - 1, ...]

            self.data.append(vals)

    def get_var_names(self):
        return [self.sanitize_string(v.tostring()) for v in
                self.dataset["name_elem_var"][:]]

    def get_nod_names(self):
        return [self.sanitize_string(v.tostring()) for v in
                self.dataset["name_nod_var"][:]]

    def sanitize_string(self, s):
        s = "".join(_ for _ in takewhile(lambda a: a in string.printable, s))
        return s
