from collections import OrderedDict

from yt.frontends.tipsy.api import TipsyDataset
from yt.testing import ParticleSelectionComparison, requires_file
from yt.utilities.answer_testing.framework import (
    data_dir_load,
    nbody_answer,
    requires_ds,
    sph_answer,
)

_fields = OrderedDict(
    [
        (("all", "particle_mass"), None),
        (("all", "particle_ones"), None),
        (("all", "particle_velocity_x"), ("all", "particle_mass")),
        (("all", "particle_velocity_y"), ("all", "particle_mass")),
        (("all", "particle_velocity_z"), ("all", "particle_mass")),
    ]
)

pkdgrav = "halo1e11_run1.00400/halo1e11_run1.00400"
pkdgrav_cosmology_parameters = dict(
    current_redshift=0.0, omega_lambda=0.728, omega_matter=0.272, hubble_constant=0.702
)
pkdgrav_kwargs = dict(
    field_dtypes={"Coordinates": "d"},
    cosmology_parameters=pkdgrav_cosmology_parameters,
    unit_base={"length": (60.0, "Mpccm/h")},
)


@requires_ds(pkdgrav, big_data=True, file_check=True)
def test_pkdgrav():
    ds = data_dir_load(pkdgrav, TipsyDataset, (), kwargs=pkdgrav_kwargs)
    yield from nbody_answer(ds, "halo1e11_run1.00400", 26847360, _fields)
    psc = ParticleSelectionComparison(ds)
    psc.run_defaults()


gasoline_dmonly = "agora_1e11.00400/agora_1e11.00400"


@requires_ds(gasoline_dmonly, big_data=True, file_check=True)
def test_gasoline_dmonly():
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
    yield from nbody_answer(ds, "agora_1e11.00400", 10550576, _fields)
    psc = ParticleSelectionComparison(ds)
    psc.run_defaults()


tg_sph_fields = OrderedDict(
    [
        (("gas", "density"), None),
        (("gas", "temperature"), None),
        (("gas", "temperature"), ("gas", "density")),
        (("gas", "velocity_magnitude"), None),
        (("gas", "Fe_fraction"), None),
    ]
)

tg_nbody_fields = OrderedDict(
    [
        (("Stars", "Metals"), None),
    ]
)

tipsy_gal = "TipsyGalaxy/galaxy.00300"


@requires_ds(tipsy_gal)
def test_tipsy_galaxy():
    ds = data_dir_load(
        tipsy_gal,
        kwargs={"bounding_box": [[-2000, 2000], [-2000, 2000], [-2000, 2000]]},
    )
    # These tests should be re-enabled.  But the holdup is that the region
    # selector does not offset by domain_left_edge, and we have inelegant
    # selection using bboxes.
    # psc = ParticleSelectionComparison(ds)
    # psc.run_defaults()
    for test in sph_answer(ds, "galaxy.00300", 315372, tg_sph_fields):
        test_tipsy_galaxy.__name__ = test.description
        yield test
    for test in nbody_answer(ds, "galaxy.00300", 315372, tg_nbody_fields):
        test_tipsy_galaxy.__name__ = test.description
        yield test


@requires_file(gasoline_dmonly)
@requires_file(pkdgrav)
def test_TipsyDataset():
    assert isinstance(data_dir_load(pkdgrav, kwargs=pkdgrav_kwargs), TipsyDataset)
    assert isinstance(data_dir_load(gasoline_dmonly), TipsyDataset)


@requires_file(tipsy_gal)
def test_tipsy_index():
    ds = data_dir_load(tipsy_gal)
    sl = ds.slice("z", 0.0)
    assert sl["gas", "density"].shape[0] != 0
