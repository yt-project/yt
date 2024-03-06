import numpy as np
import pytest
import unyt

import yt
from yt.geometry.coordinates.spherical_coordinates import spherical_to_cartesian
from yt.testing import fake_amr_ds


def test_cutting_plane_mixed_coords():
    ds = fake_amr_ds(geometry="spherical")
    normal = np.array([0.0, 0.0, 1.0])
    plane_center = np.array([0.0, 0.0, 0.5])
    slc = ds.cutting_mixed(normal, plane_center)
    frb = slc.to_frb(2.0, 800)
    bvals = frb[("index", "r")]
    mask = frb.get_mask(("index", "r"))
    # note: the min value of r on the plane will be the z value of the
    # plane center. how close it is to the correct answer will depend
    # on the size of the elements.
    assert np.allclose(bvals[mask].min().d, plane_center[2], atol=0.02)


@pytest.fixture
def spherical_ds():

    shp = (32, 32, 32)
    data = {"density": np.random.random(shp)}

    bbox = np.array([[0.0, 1.0], [0, np.pi], [0, 2 * np.pi]])

    def _z(field, data):
        r = data["index", "r"]
        theta = data["index", "theta"]
        phi = data["index", "phi"]
        _, _, z = spherical_to_cartesian(r, theta, phi)
        return unyt.unyt_array(z, r.units)

    ds = yt.load_uniform_grid(
        data,
        shp,
        bbox=bbox,
        geometry="spherical",
        axis_order=("r", "theta", "phi"),
        length_unit="m",
    )

    ds.add_field(
        name=("index", "z_val"), function=_z, sampling_type="cell", take_log=False
    )

    return ds


def test_cutting_plane_mixed_fixed_z(spherical_ds):
    ds = spherical_ds
    normal = np.array([0.0, 0.0, 1.0])
    center = np.array([0.0, 0.0, 0.5])
    slc = ds.cutting_mixed(normal, center)
    zvals = slc["index", "z_val"].to("code_length")
    assert np.allclose(zvals, ds.quan(0.5, "code_length"), atol=0.05)


@pytest.mark.answer_test
class TestCuttingPlaneMixerdCoords:
    answer_file = None
    saved_hashes = None
    answer_version = "000"

    @pytest.mark.usefixtures("hashing")
    def test_vertical_slice_at_sphere_edge(self, spherical_ds):
        ds = spherical_ds
        normal = np.array([0.0, 1.0, 0.0])
        center = np.array([0.0, 0.9, 0.0])
        slc = ds.cutting_mixed(normal, center)
        frb = slc.to_frb(2.0, 50)
        vals = frb["index", "z_val"].to("code_length")
        vals[~frb.get_mask(("index", "z_val"))] = np.nan
        vals = vals.to_ndarray()
        self.hashes.update({"offset_slice": vals})
