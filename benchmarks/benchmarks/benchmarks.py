# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
import numpy as np
from yt import YTArray, YTQuantity

def time_quantity_init_scalar1():
    3.0 * YTQuantity(1, "m/s")


def time_quantity_init_scalar2():
    YTQuantity(3.0, "m/s")


def time_quantity_init_array():
    np.arange(100000) * YTQuantity(1, "m/s")


def time_quantity_init_array2():
    YTArray(np.arange(100000), "m/s")


def time_quantity_scalar_conversion():
    YTQuantity(3.0, "m/s").in_units("km/hr")


def time_quantity_array_conversion():
    YTArray(np.arange(100000), "m/s").in_units("km/hr")


def time_quantity_ufunc_sin():
    np.sin(YTArray(np.arange(10000), "degree"))
