from yt.testing import *
import numpy as np
from yt.utilities.lib.ragged_arrays import index_unop

operations = ((np.sum, "sum"),
              (np.prod, "prod"),
              (np.max, "max"),
              (np.min, "min"))
dtypes = ((-1e8, 1e8, "float32"),
          (-1e8, 1e8, "float64"),
          (-10000, 10000, "int32"),
          (-100000000, 100000000, "int64"))

def test_index_unop():
    np.random.seed(0x4d3d3d3)
    indices = np.arange(1000, dtype="int64")
    np.random.shuffle(indices)
    sizes = np.array([
        200, 50, 50, 100, 32, 32, 32, 32, 32, 64, 376], dtype="int64")
    for mi, ma, dtype in dtypes:
        for op, operation in operations:
            # Create a random set of values
            values = np.random.random(1000)
            if operation != "prod":
                values = values * ma + (ma - mi)
            if operation == "prod" and dtype.startswith("int"):
                values = values.astype(dtype)
                values[values != 0] = 1
                values[values == 0] = -1
            values = values.astype(dtype)
            out_values = index_unop(values, indices, sizes, operation)
            i = 0
            for j, v in enumerate(sizes):
                arr = values[indices[i:i+v]]
                yield assert_equal, op(arr), out_values[j]
                i += v
