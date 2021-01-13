import pytest

from yt.frontends.tipsy.api import TipsyDataset
from yt.testing import ParticleSelectionComparison, requires_file
from yt.utilities.answer_testing.answer_tests import nbody_answer, sph_answer
from yt.utilities.answer_testing.testing_utilities import data_dir_load, requires_ds

# Test data
pkdgrav = "halo1e11_run1.00400/halo1e11_run1.00400"
gasoline_dmonly = "agora_1e11.00400/agora_1e11.00400"
tipsy_gal = "TipsyGalaxy/galaxy.00300"


pkdgrav_cosmology_parameters = dict(
    current_redshift=0.0, omega_lambda=0.728, omega_matter=0.272, hubble_constant=0.702
)

pkdgrav_kwargs = dict(
    field_dtypes={"Coordinates": "d"},
    cosmology_parameters=pkdgrav_cosmology_parameters,
    unit_base={"length": (60.0, "Mpccm/h")},
)


fields = {
    ("all", "particle_mass"): None,
    ("all", "particle_ones"): None,
    ("all", "particle_velocity_x"): ("all", "particle_mass"),
    ("all", "particle_velocity_y"): ("all", "particle_mass"),
    ("all", "particle_velocity_z"): ("all", "particle_mass"),
}

# Is a list instead of a dictionary because ("gas", "temperature")
# appeared as a duplicate key
tg_sph_fields = [
    [("gas", "density"), None],
    [("gas", "temperature"), None],
    [("gas", "temperature"), ("gas", "density")],
    [("gas", "velocity_magnitude"), None],
    [("gas", "Fe_fraction"), None],
]

tg_nbody_fields = {
    ("Stars", "Metals"): None,
}

axes = [0, 1, 2]


pkdgrav_pairs = []
ds = [
    pkdgrav,
    {
        "cls" : TipsyDataset,
        "args" : (),
        "kwargs" : pkdgrav_kwargs,
    }
]
d = [None, ("sphere", ("c", (0.3, "unitary")))]
for f, w in fields.items():
    for obj in d:
        pkdgrav_pairs.append((ds, f, obj, w))


gasoline_pairs = []
cosmology_parameters = dict(
    current_redshift=0.0,
    omega_lambda=0.728,
    omega_matter=0.272,
    hubble_constant=0.702,
)
gasoline_kwargs = dict(
    cosmology_parameters=cosmology_parameters,
    unit_base={"length": (60.0, "Mpccm/h")},
)
ds = [
    gasoline_dmonly,
    {
        "cls" : TipsyDataset,
        "args" : (),
        "kwargs" : gasoline_kwargs,
    }
]
d = [None, ("sphere", ("c", (0.3, "unitary")))]
for f, w in fields.items():
    for obj in d:
        gasoline_pairs.append((ds, f, obj, w))


galaxy_sph_pairs = []
ds = [
    tipsy_gal,
    {
        "kwargs" : {"bounding_box": [[-2000, 2000], [-2000, 2000], [-2000, 2000]]},
    },
]
d = [None, ("sphere", ("c", (0.1, "unitary")))]
for f_w in tg_sph_fields:
    for obj in d:
        galaxy_sph_pairs.append((ds, f_w[0], obj, f_w[1]))


galaxy_nbody_pairs = []
ds = [
    tipsy_gal,
    {
        "kwargs" : {"bounding_box": [[-2000, 2000], [-2000, 2000], [-2000, 2000]]},
    },
]
d = [None, ("sphere", ("c", (0.1, "unitary")))]
for f, w in tg_nbody_fields.items():
    for obj in d:
        galaxy_nbody_pairs.append((ds, f, obj, w))


@pytest.mark.answer_test
class TestTipsy:
    answer_file = None
    saved_hashes = None

    @pytest.mark.big_data
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f, d, w", pkdgrav_pairs, indirect=True)
    @pytest.mark.parametrize("a", axes, indirect=True)
    def test_pkdgrav(self, f, a, d, w, ds):
        nb = nbody_answer(ds, "halo1e11_run1.00400", 26847360, f, w, d, a)
        self.hashes.update(nb)

    @pytest.mark.big_data
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f, d, w", gasoline_pairs, indirect=True)
    @pytest.mark.parametrize("a", axes, indirect=True)
    def test_gasoline_dmonly(self, f, a, d, w, ds):
        nb = nbody_answer(ds, "agora_1e11.00400", 10550576, f, w, d, a)
        self.hashes.update(nb)

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f, d, w", galaxy_sph_pairs, indirect=True)
    @pytest.mark.parametrize("a", axes, indirect=True)
    def test_tipsy_galaxy_sph(self, f, w, d, a, ds):
        # These tests should be re-enabled.  But the holdup is that the region
        # selector does not offset by domain_left_edge, and we have inelegant
        # selection using bboxes.
        # psc = ParticleSelectionComparison(ds)
        # psc.run_defaults()
        sph = sph_answer(ds, "galaxy.00300", 315372, f, w, d, a)
        self.hashes.update(sph)

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f, d, w", galaxy_nbody_pairs, indirect=True)
    @pytest.mark.parametrize("a", axes, indirect=True)
    def test_tipsy_galaxy_nbody(self, f, w, d, a, ds):
        nb = nbody_answer(ds, "galaxy.00300", 315372, f, w, d, a)
        self.hashes.update(nb)

    @pytest.mark.big_data
    @requires_ds(pkdgrav, file_check=True)
    def test_pkdgrav_psc(self):
        ds = data_dir_load(pkdgrav, TipsyDataset, (), kwargs=pkdgrav_kwargs)
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()

    @pytest.mark.big_data
    @requires_ds(gasoline_dmonly, file_check=True)
    def test_gasoline_dmonly_psc(self):
        cosmology_parameters = dict(
            current_redshift=0.0,
            omega_lambda=0.728,
            omega_matter=0.272,
            hubble_constant=0.702,
        )
        kwargs = dict(
            cosmology_parameters=cosmology_parameters,
            unit_base={"length": (60.0, "Mpccm/h")},
        )
        ds = data_dir_load(gasoline_dmonly, TipsyDataset, (), kwargs)
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()

    @requires_file(gasoline_dmonly)
    @requires_file(pkdgrav)
    def test_TipsyDataset(self):
        assert isinstance(data_dir_load(pkdgrav, kwargs=pkdgrav_kwargs), TipsyDataset)
        assert isinstance(data_dir_load(gasoline_dmonly), TipsyDataset)

    @requires_file(tipsy_gal)
    def test_tipsy_index(self):
        ds = data_dir_load(tipsy_gal)
        sl = ds.slice("z", 0.0)
        assert sl["gas", "density"].shape[0] != 0
