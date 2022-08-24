import numpy as np

from yt import load
from yt.frontends.stream.fields import StreamFieldInfo
from yt.testing import (
    assert_allclose_units,
    assert_almost_equal,
    assert_array_almost_equal_nulp,
    assert_array_equal,
    assert_equal,
    assert_raises,
    fake_amr_ds,
    fake_particle_ds,
    fake_random_ds,
    requires_file,
)
from yt.units.yt_array import YTArray, YTQuantity, array_like_field
from yt.utilities.cosmology import Cosmology
from yt.utilities.exceptions import (
    YTDimensionalityError,
    YTFieldUnitError,
    YTFieldUnitParseError,
)


def get_params(ds):
    return dict(
        axis=0,
        center=YTArray((0.0, 0.0, 0.0), "cm", registry=ds.unit_registry),
        bulk_velocity=YTArray((0.0, 0.0, 0.0), "cm/s", registry=ds.unit_registry),
        bulk_magnetic_field=YTArray((0.0, 0.0, 0.0), "G", registry=ds.unit_registry),
        normal=YTArray((0.0, 0.0, 1.0), "", registry=ds.unit_registry),
        cp_x_vec=YTArray((1.0, 0.0, 0.0), "", registry=ds.unit_registry),
        cp_y_vec=YTArray((0.0, 1.0, 0.0), "", registry=ds.unit_registry),
        cp_z_vec=YTArray((0.0, 0.0, 1.0), "", registry=ds.unit_registry),
        omega_baryon=0.04,
        observer_redshift=0.0,
        source_redshift=3.0,
        virial_radius=YTQuantity(1.0, "cm"),
    )


_base_fields = (
    ("gas", "density"),
    ("gas", "velocity_x"),
    ("gas", "velocity_y"),
    ("gas", "velocity_z"),
)


def _strip_ftype(field):
    if not isinstance(field, tuple):
        return field
    elif field[0] in ("all", "io"):
        return field
    return field[1]


class TestFieldAccess:
    description = None

    def __init__(self, field_name, ds, nprocs):
        # Note this should be a field name
        self.field_name = field_name
        self.description = f"Accessing_{field_name}_{nprocs}"
        self.nprocs = nprocs
        self.ds = ds

    def __call__(self):
        field = self.ds._get_field_info(*self.field_name)
        skip_grids = False
        needs_spatial = False
        for v in field.validators:
            if getattr(v, "ghost_zones", 0) > 0:
                skip_grids = True
            if hasattr(v, "ghost_zones"):
                needs_spatial = True

        ds = self.ds

        # This gives unequal sized grids as well as subgrids
        dd1 = ds.all_data()
        dd2 = ds.all_data()
        sp = get_params(ds)
        dd1.field_parameters.update(sp)
        dd2.field_parameters.update(sp)
        with np.errstate(all="ignore"):
            v1 = dd1[self.field_name]
            # No more conversion checking
            assert_equal(v1, dd1[self.field_name])
            if not needs_spatial:
                with field.unit_registry(dd2):
                    res = field._function(field, dd2)
                    res = dd2.apply_units(res, field.units)
                assert_array_almost_equal_nulp(v1, res, 4)
            if not skip_grids:
                for g in ds.index.grids:
                    g.field_parameters.update(sp)
                    v1 = g[self.field_name]
                    g.clear_data()
                    g.field_parameters.update(sp)
                    r1 = field._function(field, g)
                    if field.sampling_type == "particle":
                        assert_equal(v1.shape[0], g.NumberOfParticles)
                    else:
                        assert_array_equal(r1.shape, v1.shape)
                        for ax in "xyz":
                            assert_array_equal(g["index", ax].shape, v1.shape)
                    with field.unit_registry(g):
                        res = field._function(field, g)
                        assert_array_equal(v1.shape, res.shape)
                        res = g.apply_units(res, field.units)
                    assert_array_almost_equal_nulp(v1, res, 4)


def get_base_ds(nprocs):
    fields, units = [], []

    for fname, (code_units, *_) in StreamFieldInfo.known_other_fields:
        fields.append(("gas", fname))
        units.append(code_units)

    pfields, punits = [], []

    for fname, (code_units, _aliases, _dn) in StreamFieldInfo.known_particle_fields:
        if fname == "smoothing_lenth":
            # we test SPH fields elsewhere
            continue
        pfields.append(fname)
        punits.append(code_units)

    ds = fake_random_ds(
        4,
        fields=fields,
        units=units,
        particles=20,
        nprocs=nprocs,
        particle_fields=pfields,
        particle_field_units=punits,
    )
    ds.parameters["HydroMethod"] = "streaming"
    ds.parameters["EOSType"] = 1.0
    ds.parameters["EOSSoundSpeed"] = 1.0
    ds.conversion_factors["Time"] = 1.0
    ds.conversion_factors.update({f: 1.0 for f in fields})
    ds.gamma = 5.0 / 3.0
    ds.current_redshift = 0.0001
    ds.cosmological_simulation = 1
    ds.hubble_constant = 0.7
    ds.omega_matter = 0.27
    ds.omega_lambda = 0.73
    ds.cosmology = Cosmology(
        hubble_constant=ds.hubble_constant,
        omega_matter=ds.omega_matter,
        omega_lambda=ds.omega_lambda,
        unit_registry=ds.unit_registry,
    )
    # ensures field errors are raised during testing
    # see FieldInfoContainer.check_derived_fields
    ds._field_test_dataset = True
    ds.index
    return ds


def test_all_fields():
    datasets = {}

    for nprocs in [1, 4, 8]:
        ds = get_base_ds(nprocs)
        datasets[nprocs] = ds

    for field in sorted(ds.field_info):
        if field[1].find("beta_p") > -1:
            continue
        if field[1].find("vertex") > -1:
            # don't test the vertex fields for now
            continue
        if field[1].find("smoothed") > -1:
            # smoothed fields aren't implemented for grid data
            continue
        if field in ds.field_list:
            # Don't know how to test this.  We need some way of having fields
            # that are fallbacks be tested, but we don't have that now.
            continue

        for nprocs in [1, 4, 8]:
            test_all_fields.__name__ = f"{field}_{nprocs}"
            yield TestFieldAccess(field, datasets[nprocs], nprocs)


def test_add_deposited_particle_field():
    # NOT tested: "std", "mesh_id", "nearest" and "simple_smooth"
    base_ds = get_base_ds(1)
    ad = base_ds.all_data()

    # Test "count", "sum" and "cic" method
    for method in ["count", "sum", "cic"]:
        fn = base_ds.add_deposited_particle_field(("io", "particle_mass"), method)
        expected_fn = "io_%s" if method == "count" else "io_%s_mass"
        assert_equal(fn, ("deposit", expected_fn % method))
        ret = ad[fn]
        if method == "count":
            assert_equal(ret.sum(), ad[("io", "particle_ones")].sum())
        else:
            assert_almost_equal(ret.sum(), ad[("io", "particle_mass")].sum())

    # Test "weighted_mean" method
    fn = base_ds.add_deposited_particle_field(
        ("io", "particle_ones"), "weighted_mean", weight_field="particle_ones"
    )
    assert_equal(fn, ("deposit", "io_avg_ones"))
    ret = ad[fn]
    # The sum should equal the number of cells that have particles
    assert_equal(ret.sum(), np.count_nonzero(ad[("deposit", "io_count")]))


def test_add_gradient_fields():
    ds = get_base_ds(1)
    gfields = ds.add_gradient_fields(("gas", "density"))
    gfields += ds.add_gradient_fields(("index", "ones"))
    field_list = [
        ("gas", "density_gradient_x"),
        ("gas", "density_gradient_y"),
        ("gas", "density_gradient_z"),
        ("gas", "density_gradient_magnitude"),
        ("index", "ones_gradient_x"),
        ("index", "ones_gradient_y"),
        ("index", "ones_gradient_z"),
        ("index", "ones_gradient_magnitude"),
    ]
    assert_equal(gfields, field_list)
    ad = ds.all_data()
    for field in field_list:
        ret = ad[field]
        if field[0] == "gas":
            assert str(ret.units) == "g/cm**4"
        else:
            assert str(ret.units) == "1/cm"


def test_add_gradient_fields_by_fname():
    ds = fake_amr_ds(fields=("density", "temperature"), units=("g/cm**3", "K"))
    actual = ds.add_gradient_fields(("gas", "density"))
    expected = [
        ("gas", "density_gradient_x"),
        ("gas", "density_gradient_y"),
        ("gas", "density_gradient_z"),
        ("gas", "density_gradient_magnitude"),
    ]
    assert_equal(actual, expected)


def test_add_gradient_multiple_fields():
    ds = fake_amr_ds(fields=("density", "temperature"), units=("g/cm**3", "K"))
    actual = ds.add_gradient_fields([("gas", "density"), ("gas", "temperature")])
    expected = [
        ("gas", "density_gradient_x"),
        ("gas", "density_gradient_y"),
        ("gas", "density_gradient_z"),
        ("gas", "density_gradient_magnitude"),
        ("gas", "temperature_gradient_x"),
        ("gas", "temperature_gradient_y"),
        ("gas", "temperature_gradient_z"),
        ("gas", "temperature_gradient_magnitude"),
    ]
    assert_equal(actual, expected)

    ds = fake_amr_ds(fields=("density", "temperature"), units=("g/cm**3", "K"))
    actual = ds.add_gradient_fields([("gas", "density"), ("gas", "temperature")])
    assert_equal(actual, expected)


def test_add_gradient_fields_curvilinear():
    ds = fake_amr_ds(fields=["density"], units=["g/cm**3"], geometry="spherical")
    gfields = ds.add_gradient_fields(("gas", "density"))
    gfields += ds.add_gradient_fields(("index", "ones"))
    field_list = [
        ("gas", "density_gradient_r"),
        ("gas", "density_gradient_theta"),
        ("gas", "density_gradient_phi"),
        ("gas", "density_gradient_magnitude"),
        ("index", "ones_gradient_r"),
        ("index", "ones_gradient_theta"),
        ("index", "ones_gradient_phi"),
        ("index", "ones_gradient_magnitude"),
    ]
    assert_equal(gfields, field_list)
    ad = ds.all_data()
    for field in field_list:
        ret = ad[field]
        if field[0] == "gas":
            assert str(ret.units) == "g/cm**4"
        else:
            assert str(ret.units) == "1/cm"


def get_data(ds, field_name):
    # Need to create a new data object otherwise the errors we are
    # intentionally raising lead to spurious GenerationInProgress errors
    ad = ds.all_data()
    return ad[field_name]


def test_add_field_unit_semantics():
    ds = fake_random_ds(16)
    ad = ds.all_data()

    def density_alias(field, data):
        return data[("gas", "density")].in_cgs()

    def unitless_data(field, data):
        return np.ones(data[("gas", "density")].shape)

    ds.add_field(
        ("gas", "density_alias_auto"),
        sampling_type="cell",
        function=density_alias,
        units="auto",
        dimensions="density",
    )
    ds.add_field(
        ("gas", "density_alias_wrong_units"),
        function=density_alias,
        sampling_type="cell",
        units="m/s",
    )
    ds.add_field(
        ("gas", "density_alias_unparseable_units"),
        sampling_type="cell",
        function=density_alias,
        units="dragons",
    )
    ds.add_field(
        ("gas", "density_alias_auto_wrong_dims"),
        function=density_alias,
        sampling_type="cell",
        units="auto",
        dimensions="temperature",
    )
    assert_raises(YTFieldUnitError, get_data, ds, ("gas", "density_alias_wrong_units"))
    assert_raises(
        YTFieldUnitParseError, get_data, ds, ("gas", "density_alias_unparseable_units")
    )
    assert_raises(
        YTDimensionalityError, get_data, ds, ("gas", "density_alias_auto_wrong_dims")
    )

    dens = ad[("gas", "density_alias_auto")]
    assert_equal(str(dens.units), "g/cm**3")

    ds.add_field(("gas", "dimensionless"), sampling_type="cell", function=unitless_data)
    ds.add_field(
        ("gas", "dimensionless_auto"),
        function=unitless_data,
        sampling_type="cell",
        units="auto",
        dimensions="dimensionless",
    )
    ds.add_field(
        ("gas", "dimensionless_explicit"),
        function=unitless_data,
        sampling_type="cell",
        units="",
    )
    ds.add_field(
        ("gas", "dimensionful"),
        sampling_type="cell",
        function=unitless_data,
        units="g/cm**3",
    )

    assert_equal(str(ad[("gas", "dimensionless")].units), "dimensionless")
    assert_equal(str(ad[("gas", "dimensionless_auto")].units), "dimensionless")
    assert_equal(str(ad[("gas", "dimensionless_explicit")].units), "dimensionless")
    assert_raises(YTFieldUnitError, get_data, ds, ("gas", "dimensionful"))


def test_add_field_from_lambda():
    ds = fake_amr_ds(fields=["density"], units=["g/cm**3"])

    def _function_density(field, data):
        return data["gas", "density"]

    ds.add_field(
        ("gas", "function_density"),
        function=_function_density,
        sampling_type="cell",
        units="g/cm**3",
    )

    ds.add_field(
        ("gas", "lambda_density"),
        function=lambda field, data: data["gas", "density"],
        sampling_type="cell",
        units="g/cm**3",
    )

    ad = ds.all_data()
    # check that the fields are accessible
    ad["gas", "function_density"]
    ad["gas", "lambda_density"]


def test_array_like_field():
    ds = fake_random_ds(4, particles=64)
    ad = ds.all_data()
    u1 = ad[("all", "particle_mass")].units
    u2 = array_like_field(ad, 1.0, ("all", "particle_mass")).units
    assert u1 == u2


ISOGAL = "IsolatedGalaxy/galaxy0030/galaxy0030"


@requires_file(ISOGAL)
def test_array_like_field_output_units():
    ds = load(ISOGAL)
    ad = ds.all_data()
    u1 = ad[("all", "particle_mass")].units
    u2 = array_like_field(ad, 1.0, ("all", "particle_mass")).units
    assert u1 == u2
    assert str(u1) == ds.fields.all.particle_mass.output_units
    u1 = ad[("gas", "x")].units
    u2 = array_like_field(ad, 1.0, ("gas", "x")).units
    assert u1 == u2
    assert str(u1) == ds.fields.gas.x.units


def test_add_field_string():
    ds = fake_random_ds(16)
    ad = ds.all_data()

    def density_alias(field, data):
        return data[("gas", "density")]

    ds.add_field(
        ("gas", "density_alias"),
        sampling_type="cell",
        function=density_alias,
        units="g/cm**3",
    )

    ad[("gas", "density_alias")]

    assert ("gas", "density_alias") in ds.derived_field_list


def test_add_field_string_aliasing():
    ds = fake_random_ds(16)

    def density_alias(field, data):
        return data["gas", "density"]

    ds.add_field(
        ("gas", "density_alias"),
        sampling_type="cell",
        function=density_alias,
        units="g/cm**3",
    )

    ds.field_info["gas", "density_alias"]
    ds.field_info["gas", "density_alias"]

    ds = fake_particle_ds()

    def pmass_alias(field, data):
        return data["all", "particle_mass"]

    ds.add_field(
        ("all", "particle_mass_alias"),
        function=pmass_alias,
        units="g",
        sampling_type="particle",
    )

    ds.field_info["all", "particle_mass_alias"]
    ds.field_info["all", "particle_mass_alias"]


def test_morton_index():
    ds = fake_amr_ds()
    mi = ds.r["index", "morton_index"]
    mi2 = mi.view("uint64")
    assert_equal(np.unique(mi2).size, mi2.size)
    a1 = np.argsort(mi)
    a2 = np.argsort(mi2)
    assert_array_equal(a1, a2)


def test_field_inference():
    ds = fake_random_ds(16)
    ds.index
    # If this is not true this means the result of field inference depends
    # on the order we did field detection, which is random in Python3
    assert_equal(ds._last_freq, (None, None))


@requires_file(ISOGAL)
def test_deposit_amr():
    ds = load(ISOGAL)
    for g in ds.index.grids:
        gpm = g["all", "particle_mass"].sum()
        dpm = g["deposit", "all_mass"].sum()
        assert_allclose_units(gpm, dpm)


def test_ion_field_labels():
    fields = [
        "O_p1_number_density",
        "O2_p1_number_density",
        "CO2_p1_number_density",
        "Co_p1_number_density",
        "O2_p2_number_density",
        "H2O_p1_number_density",
    ]
    units = ["cm**-3" for f in fields]
    ds = fake_random_ds(16, fields=fields, units=units)

    # by default labels should use roman numerals
    default_labels = {
        "O_p1_number_density": "$\\rm{O\\ II\\ Number\\ Density}$",
        "O2_p1_number_density": "$\\rm{O_{2}\\ II\\ Number\\ Density}$",
        "CO2_p1_number_density": "$\\rm{CO_{2}\\ II\\ Number\\ Density}$",
        "Co_p1_number_density": "$\\rm{Co\\ II\\ Number\\ Density}$",
        "O2_p2_number_density": "$\\rm{O_{2}\\ III\\ Number\\ Density}$",
        "H2O_p1_number_density": "$\\rm{H_{2}O\\ II\\ Number\\ Density}$",
    }

    pm_labels = {
        "O_p1_number_density": "$\\rm{{O}^{+}\\ Number\\ Density}$",
        "O2_p1_number_density": "$\\rm{{O_{2}}^{+}\\ Number\\ Density}$",
        "CO2_p1_number_density": "$\\rm{{CO_{2}}^{+}\\ Number\\ Density}$",
        "Co_p1_number_density": "$\\rm{{Co}^{+}\\ Number\\ Density}$",
        "O2_p2_number_density": "$\\rm{{O_{2}}^{++}\\ Number\\ Density}$",
        "H2O_p1_number_density": "$\\rm{{H_{2}O}^{+}\\ Number\\ Density}$",
    }

    fobj = ds.fields.stream

    for f in fields:
        label = getattr(fobj, f).get_latex_display_name()
        assert_equal(label, default_labels[f])

    ds.set_field_label_format("ionization_label", "plus_minus")
    fobj = ds.fields.stream

    for f in fields:
        label = getattr(fobj, f).get_latex_display_name()
        assert_equal(label, pm_labels[f])


def test_default_fluid_type_None():
    """
    Check for bad behavior when default_fluid_type is None.
    See PR #3710.
    """
    ds = fake_amr_ds()
    ds.default_fluid_type = None
    ds.field_list
