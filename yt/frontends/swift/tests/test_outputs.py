import numpy as np
from numpy.testing import assert_almost_equal

from yt import load
from yt.frontends.swift.api import SwiftDataset
from yt.testing import ParticleSelectionComparison, requires_file, requires_module
from yt.utilities.on_demand_imports import _h5py as h5py

keplerian_ring = "KeplerianRing/keplerian_ring_0020.hdf5"
EAGLE_6 = "EAGLE_6/eagle_0005.hdf5"
stromgen_sphere = "Stromgren/output_7bin_HHeZdust_real_0060.hdf5"


# Combined the tests for loading a file and ensuring the units have been
# implemented correctly to save time on re-loading a dataset
@requires_module("h5py")
@requires_file(keplerian_ring)
def test_non_cosmo_dataset():
    ds = load(keplerian_ring)
    assert type(ds) is SwiftDataset

    field = ("gas", "density")
    ad = ds.all_data()
    yt_density = ad[field]
    yt_coords = ad[field[0], "position"]

    # load some data the old fashioned way
    fh = h5py.File(ds.parameter_filename, mode="r")
    part_data = fh["PartType0"]

    # set up a conversion factor by loading the unit mas and unit length in cm,
    # and then converting to proper coordinates
    units = fh["Units"]
    units = dict(units.attrs)
    density_factor = float(units["Unit mass in cgs (U_M)"])
    density_factor /= float(units["Unit length in cgs (U_L)"]) ** 3

    # now load the raw density and coordinates
    raw_density = part_data["Density"][:].astype("float64") * density_factor
    raw_coords = part_data["Coordinates"][:].astype("float64")
    fh.close()

    # sort by the positions - yt often loads in a different order
    ind_raw = np.lexsort((raw_coords[:, 2], raw_coords[:, 1], raw_coords[:, 0]))
    ind_yt = np.lexsort((yt_coords[:, 2], yt_coords[:, 1], yt_coords[:, 0]))
    raw_density = raw_density[ind_raw]
    yt_density = yt_density[ind_yt]

    # make sure we are comparing fair units
    assert str(yt_density.units) == "g/cm**3"

    # make sure the actual values are the same
    assert_almost_equal(yt_density.d, raw_density)


@requires_module("h5py")
@requires_file(keplerian_ring)
def test_non_cosmo_dataset_selection():
    ds = load(keplerian_ring)
    psc = ParticleSelectionComparison(ds)
    psc.run_defaults()


@requires_module("h5py")
@requires_file(EAGLE_6)
def test_cosmo_dataset():
    ds = load(EAGLE_6)
    assert type(ds) is SwiftDataset

    field = ("gas", "density")
    ad = ds.all_data()
    yt_density = ad[field]
    yt_coords = ad[field[0], "position"]

    # load some data the old fashioned way
    fh = h5py.File(ds.parameter_filename, mode="r")
    part_data = fh["PartType0"]

    # set up a conversion factor by loading the unit mas and unit length in cm,
    # and then converting to proper coordinates
    units = fh["Units"]
    units = dict(units.attrs)
    density_factor = float(units["Unit mass in cgs (U_M)"])
    density_factor /= float(units["Unit length in cgs (U_L)"]) ** 3

    # add the redshift factor
    header = fh["Header"]
    header = dict(header.attrs)
    density_factor *= (1.0 + float(header["Redshift"])) ** 3

    # now load the raw density and coordinates
    raw_density = part_data["Density"][:].astype("float64") * density_factor
    raw_coords = part_data["Coordinates"][:].astype("float64")
    fh.close()

    # sort by the positions - yt often loads in a different order
    ind_raw = np.lexsort((raw_coords[:, 2], raw_coords[:, 1], raw_coords[:, 0]))
    ind_yt = np.lexsort((yt_coords[:, 2], yt_coords[:, 1], yt_coords[:, 0]))
    raw_density = raw_density[ind_raw]
    yt_density = yt_density[ind_yt]

    # make sure we are comparing fair units
    assert str(yt_density.units) == "g/cm**3"

    # make sure the actual values are the same
    assert_almost_equal(yt_density.d, raw_density)


@requires_module("h5py")
@requires_file(EAGLE_6)
def test_cosmo_dataset_selection():
    ds = load(EAGLE_6)
    psc = ParticleSelectionComparison(ds)
    psc.run_defaults()


@requires_module("h5py")
@requires_file(stromgen_sphere)
def test_species_fractions():
    """
    Chemistry-network SWIFT snapshots expose a packed "SpeciesFractions"
    field; SwiftFieldInfo unpack it into individual per-ion and
    per-molecule "<species>_fraction" (normalized) and
    "<species>_fraction_raw" (on-disk, unnormalized) fields.
    """
    ds = load(stromgen_sphere)
    if ("PartType0", "SpeciesFractions") not in ds.field_list:
        return  # this sample dataset has no chemistry network fields

    ad = ds.all_data()

    # normalized ion fraction always lies in [0, 1]
    for name in ("HI", "HII"):
        frac = ad["PartType0", f"{name}_fraction"]
        assert frac.min() >= 0.0
        assert frac.max() <= 1.0

    # the raw column value need not be bounded the same way, and for
    # H (single-ion-state indices aside) should be present and finite
    raw = ad["PartType0", "HII_fraction_raw"]
    assert np.isfinite(raw).all()

    # HI + HII + Hm should renormalize to 1 (mask out zero-total particles)
    hi = ad["PartType0", "HI_fraction"]
    hii = ad["PartType0", "HII_fraction"]
    hm = ad["PartType0", "Hm_fraction"]
    total_raw = (
        ad["PartType0", "HI_fraction_raw"]
        + ad["PartType0", "HII_fraction_raw"]
        + ad["PartType0", "Hm_fraction_raw"]
    )
    nonzero = total_raw > 0
    assert_almost_equal(
        (hi + hii + hm)[nonzero], np.ones(nonzero.sum(), dtype="float64")
    )

    # a molecule's normalized fraction is just an alias for its raw value
    assert_almost_equal(
        ad["PartType0", "H2_fraction"].d, ad["PartType0", "H2_fraction_raw"].d
    )
