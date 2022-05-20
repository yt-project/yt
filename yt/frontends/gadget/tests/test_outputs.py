import os
import shutil
import tempfile
from collections import OrderedDict
from itertools import product

import yt
from yt.frontends.gadget.api import GadgetDataset, GadgetHDF5Dataset
from yt.frontends.gadget.testing import fake_gadget_binary
from yt.testing import ParticleSelectionComparison, requires_file
from yt.utilities.answer_testing.framework import data_dir_load, requires_ds, sph_answer

isothermal_h5 = "IsothermalCollapse/snap_505.hdf5"
isothermal_bin = "IsothermalCollapse/snap_505"
BE_Gadget = "BigEndianGadgetBinary/BigEndianGadgetBinary"
LE_SnapFormat2 = "Gadget3-snap-format2/Gadget3-snap-format2"
keplerian_ring = "KeplerianRing/keplerian_ring_0020.hdf5"
snap_33 = "snapshot_033/snap_033.0.hdf5"
snap_33_dir = "snapshot_033/"
magneticum = "MagneticumCluster/snap_132"

# This maps from field names to weight field names to use for projections
iso_fields = OrderedDict(
    [
        (("gas", "density"), None),
        (("gas", "temperature"), None),
        (("gas", "temperature"), ("gas", "density")),
        (("gas", "velocity_magnitude"), None),
    ]
)
iso_kwargs = dict(bounding_box=[[-3, 3], [-3, 3], [-3, 3]])


def test_gadget_binary():
    header_specs = ["default", "default+pad32", ["default", "pad32"]]
    curdir = os.getcwd()
    tmpdir = tempfile.mkdtemp()
    for header_spec, endian, fmt in product(header_specs, "<>", [1, 2]):
        try:
            fake_snap = fake_gadget_binary(
                header_spec=header_spec, endian=endian, fmt=fmt
            )
        except FileNotFoundError:
            # sometimes this happens for mysterious reasons
            pass
        ds = yt.load(fake_snap, header_spec=header_spec)
        assert isinstance(ds, GadgetDataset)
        ds.field_list
        try:
            os.remove(fake_snap)
        except FileNotFoundError:
            # sometimes this happens for mysterious reasons
            pass
    os.chdir(curdir)
    shutil.rmtree(tmpdir)


@requires_file(isothermal_h5)
def test_gadget_hdf5():
    assert isinstance(
        data_dir_load(isothermal_h5, kwargs=iso_kwargs), GadgetHDF5Dataset
    )


@requires_file(keplerian_ring)
def test_non_cosmo_dataset():
    """
    Non-cosmological datasets may not have the cosmological parameters in the
    Header. The code should fall back gracefully when they are not present,
    with the Redshift set to 0.
    """
    data = data_dir_load(keplerian_ring)
    assert data.current_redshift == 0.0
    assert data.cosmological_simulation == 0


@requires_ds(isothermal_h5)
def test_iso_collapse():
    ds = data_dir_load(isothermal_h5, kwargs=iso_kwargs)
    for test in sph_answer(ds, "snap_505", 2**17, iso_fields):
        test_iso_collapse.__name__ = test.description
        yield test


@requires_ds(LE_SnapFormat2)
def test_pid_uniqueness():
    """
    ParticleIDs should be unique.
    """
    ds = data_dir_load(LE_SnapFormat2)
    ad = ds.all_data()
    pid = ad[("all", "ParticleIDs")]
    assert len(pid) == len(set(pid.v))


@requires_file(snap_33)
@requires_file(snap_33_dir)
def test_multifile_read():
    """
    Tests to make sure multi-file gadget snapshot can be loaded by passing '.0' file
    or by passing the directory containing the multi-file snapshot.
    """
    assert isinstance(data_dir_load(snap_33), GadgetDataset)
    assert isinstance(data_dir_load(snap_33_dir), GadgetDataset)


@requires_file(snap_33)
def test_particle_subselection():
    # This checks that we correctly subselect from a dataset, first by making
    # sure we get all the particles, then by comparing manual selections against
    # them.
    ds = data_dir_load(snap_33)
    psc = ParticleSelectionComparison(ds)
    psc.run_defaults()


@requires_ds(BE_Gadget)
def test_bigendian_field_access():
    ds = data_dir_load(BE_Gadget)
    data = ds.all_data()
    data["Halo", "Velocities"]


mag_fields = OrderedDict(
    [
        (("gas", "density"), None),
        (("gas", "temperature"), None),
        (("gas", "temperature"), ("gas", "density")),
        (("gas", "velocity_magnitude"), None),
        (("gas", "H_fraction"), None),
        (("gas", "C_fraction"), None),
    ]
)


mag_kwargs = dict(
    long_ids=True,
    field_spec="magneticum_box2_hr",
)


@requires_ds(magneticum)
def test_magneticum():
    ds = data_dir_load(magneticum, kwargs=mag_kwargs)
    for test in sph_answer(ds, "snap_132", 3718111, mag_fields, center="max"):
        test_magneticum.__name__ = test.description
        yield test
