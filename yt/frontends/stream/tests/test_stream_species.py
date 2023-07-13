import numpy as np

from yt.loaders import load_uniform_grid
from yt.testing import assert_allclose_units


def test_stream_species():
    prng = np.random.default_rng(seed=42)
    arr = prng.uniform(size=(32, 32, 32))

    data = {
        "density": (arr, "g/cm**3"),
        "H_p0_fraction": (0.37 * np.ones_like(arr), "dimensionless"),
        "H_p1_fraction": (0.37 * np.ones_like(arr), "dimensionless"),
        "He_fraction": (0.24 * np.ones_like(arr), "dimensionless"),
        "CO_fraction": (0.02 * np.ones_like(arr), "dimensionless"),
    }
    bbox = np.array([[-1.5, 1.5], [-1.5, 1.5], [-1.5, 1.5]])
    ds = load_uniform_grid(data, arr.shape, length_unit="Mpc", bbox=bbox, nprocs=64)
    assert ("gas", "CO_density") in ds.derived_field_list
    assert ("gas", "H_nuclei_density") in ds.derived_field_list
    assert ("gas", "H_p0_number_density") in ds.derived_field_list
    dd = ds.all_data()
    assert_allclose_units(dd["gas", "CO_density"], 0.02 * dd["gas", "density"])
    all_H = dd["gas", "H_p0_number_density"] + dd["gas", "H_p1_number_density"]
    assert_allclose_units(all_H, dd["gas", "H_nuclei_density"])
