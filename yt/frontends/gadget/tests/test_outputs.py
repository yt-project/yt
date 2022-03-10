import os
import shutil
import tempfile
from collections import OrderedDict
from itertools import product

from numpy import array_equiv

import yt
from yt.frontends.gadget.api import GadgetDataset, GadgetHDF5Dataset
from yt.frontends.gadget.data_structures import GadgetHDF5File
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
    for test in sph_answer(ds, "snap_505", 2 ** 17, iso_fields):
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


@requires_file(snap_33)
def test_particle_file_operations():
    # This runs some unit test operations related particle io
    ds = data_dir_load(snap_33)
    ad = ds.all_data()
    ad.index._identify_base_chunk(ad)
    chunks = ad.index._chunk_io(ad, cache=False)
    chunks = list(chunks)

    # check that the chunk emitting runs
    file_set = ds.index.io._data_files_set(chunks)
    assert all([isinstance(i, GadgetHDF5File) for i in file_set])

    # check that our ordering works
    pfiles_1 = list(ds.index.io._sorted_chunk_iterator(chunks))
    reverse_chunks = chunks[::-1]  # reversed input should still yield in order
    pfiles_2 = list(ds.index.io._sorted_chunk_iterator(reverse_chunks))
    assert all([p1 == p2 for p1, p2 in zip(pfiles_1, pfiles_2)])

    # check that the field sanitization returns expected field and column
    particle_file = list(file_set)[0]
    for fld in particle_file._fields_with_cols:
        fld_col = particle_file._sanitize_field_col(f"{fld}1")
        if "Chemistry" in fld:
            assert fld_col == ("ChemistryAbundances", 1)
        else:
            assert fld_col == (fld.rstrip("_"), 1)
    fld_col = particle_file._sanitize_field_col("FakeField_1")
    assert fld_col == ("FakeField_1", None)  # not known, will not find a col
    for fld in ds.index.io._element_names:
        fld_col = particle_file._sanitize_field_col(fld)
        assert fld_col == (f"ElementAbundance/{fld}", None)

    # check that we can read in a field with all permutations of the
    # transaction handle
    ptype, field = "PartType0", "Density"
    # with no existing handle (will open and close internally)
    den0 = particle_file._read_field(ptype, field)
    # with the transaction context manager
    with particle_file.transaction() as handle:
        den1 = particle_file._read_field(ptype, field, handle=handle)
        assert array_equiv(den0, den1)
        # a nested transaction will avoid re-opening
        with particle_file.transaction(handle) as f:
            assert handle == f
            den1 = particle_file._read_field(ptype, field, handle=f)
    assert array_equiv(den0, den1)


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
