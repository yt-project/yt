from yt.testing import *
import numpy as np
from yt.utilities.cosmology import \
     Cosmology
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion
from yt.frontends.stream.fields import \
    StreamFieldInfo
from yt.units.yt_array import \
     YTArray, YTQuantity

def setup():
    global base_ds
    # Make this super teeny tiny
    fields, units = [], []

    for fname, (code_units, aliases, dn) in StreamFieldInfo.known_other_fields:
        fields.append(("gas", fname))
        units.append(code_units)

    base_ds = fake_random_ds(4, fields=fields, units=units, particles=10)

    base_ds.index
    base_ds.cosmological_simulation = 1
    base_ds.cosmology = Cosmology()
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"
    np.seterr(all = 'ignore')

def get_params(ds):
    return dict(
        axis = 0,
        center = YTArray((0.0, 0.0, 0.0), "cm",
            registry = ds.unit_registry),
        bulk_velocity = YTArray((0.0, 0.0, 0.0),
            "cm/s", registry = ds.unit_registry),
        normal = YTArray((0.0, 0.0, 1.0),
            "", registry = ds.unit_registry),
        cp_x_vec = YTArray((1.0, 0.0, 0.0),
            "", registry = ds.unit_registry),
        cp_y_vec = YTArray((0.0, 1.0, 0.0),
            "", registry = ds.unit_registry),
        cp_z_vec = YTArray((0.0, 0.0, 1.0),
            "", registry = ds.unit_registry),
        omega_baryon = 0.04,
        observer_redshift = 0.0,
        source_redshift = 3.0,
        virial_radius = YTQuantity(1.0, "cm"),
    )

_base_fields = (("gas", "density"),
                ("gas", "velocity_x"),
                ("gas", "velocity_y"),
                ("gas", "velocity_z"))

def realistic_ds(fields, particle_fields, nprocs):
    np.random.seed(int(0x4d3d3d3))
    units = [base_ds._get_field_info(*f).units for f in fields]
    punits = [base_ds._get_field_info('io', f).units for f in particle_fields]
    fields = [_strip_ftype(f) for f in fields]

    ds = fake_random_ds(16, fields=fields, units=units, nprocs=nprocs,
                        particle_fields=particle_fields,
                        particle_field_units=punits,
                        particles=base_ds.stream_handler.particle_count[0][0])

    ds.parameters["HydroMethod"] = "streaming"
    ds.parameters["EOSType"] = 1.0
    ds.parameters["EOSSoundSpeed"] = 1.0
    ds.conversion_factors["Time"] = 1.0
    ds.conversion_factors.update( dict((f, 1.0) for f in fields) )
    ds.gamma = 5.0/3.0
    ds.current_redshift = 0.0001
    ds.cosmological_simulation = 1
    ds.hubble_constant = 0.7
    ds.omega_matter = 0.27
    ds.omega_lambda = 0.73
    ds.cosmology = Cosmology(hubble_constant=ds.hubble_constant,
                             omega_matter=ds.omega_matter,
                             omega_lambda=ds.omega_lambda,
                             unit_registry=ds.unit_registry)
    return ds

def _strip_ftype(field):
    if not isinstance(field, tuple):
        return field
    elif field[0] in ("all", "io"):
        return field
    return field[1]

def _expand_field(field):
    if isinstance(field, tuple):
        return field
    if field in KnownStreamFields:
        fi = KnownStreamFields[field]
        if fi.particle_type:
            return ("all", field)
        else:
            return ("gas", field)
    # Otherwise, we just guess.
    if "particle" in field:
        return ("all", field)
    return ("gas", field)

class TestFieldAccess(object):
    description = None

    def __init__(self, field_name, nproc):
        # Note this should be a field name
        self.field_name = field_name
        self.description = "Accessing_%s_%s" % (field_name, nproc)
        self.nproc = nproc

    def __call__(self):

        field = base_ds._get_field_info(*self.field_name)
        deps = field.get_dependencies(ds = base_ds)
        requested = deps.requested
        particle_fields = \
            ['particle_position_x', 'particle_position_y', 'particle_position_z',
             'particle_velocity_x', 'particle_velocity_y', 'particle_velocity_z',
             'particle_mass']
        fields = list(_base_fields)

        for rf in requested:
            if field.particle_type:
                if rf not in particle_fields:
                    particle_fields.append(rf[1])
            else:
                fields.append(rf)

        skip_grids = False
        needs_spatial = False
        for v in field.validators:
            f = getattr(v, "fields", None)
            if f: fields += f
            if getattr(v, "ghost_zones", 0) > 0:
                skip_grids = True
            if hasattr(v, "ghost_zones"):
                needs_spatial = True

        ds = realistic_ds(fields, particle_fields, self.nproc)

        # This gives unequal sized grids as well as subgrids
        dd1 = ds.all_data()
        dd2 = ds.all_data()
        sp = get_params(ds)
        dd1.field_parameters.update(sp)
        dd2.field_parameters.update(sp)
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
                if field.particle_type:
                    assert_equal(v1.shape[0], g.NumberOfParticles)
                else:
                    assert_array_equal(r1.shape, v1.shape)
                    for ax in 'xyz':
                        assert_array_equal(g[ax].shape, v1.shape)
                with field.unit_registry(g):
                    res = field._function(field, g)
                    assert_array_equal(v1.shape, res.shape)
                    res = g.apply_units(res, field.units)
                assert_array_almost_equal_nulp(v1, res, 4)

def test_all_fields():
    for field in sorted(base_ds.field_info):
        if field[1].find("beta_p") > -1:
            continue
        if field in base_ds.field_list:
            # Don't know how to test this.  We need some way of having fields
            # that are fallbacks be tested, but we don't have that now.
            continue

        for nproc in [1, 4, 8]:
            test_all_fields.__name__ = "%s_%s" % (field, nproc)
            yield TestFieldAccess(field, nproc)

def test_add_deposited_particle_field():
    fn = base_ds.add_deposited_particle_field(('io', 'particle_ones'), 'count')
    assert_equal(fn, ('deposit', 'io_count_ones'))
    ad = base_ds.all_data()
    ret = ad[fn]
    assert_equal(ret.sum(), ad['particle_ones'].sum())

def get_data(ds, field_name):
    # Need to create a new data object otherwise the errors we are
    # intentionally raising lead to spurious GenerationInProgress errors
    ad = ds.all_data()
    return ad[field_name]

def test_add_field_unit_semantics():
    ds = fake_random_ds(16)
    ad = ds.all_data()

    def density_alias(field, data):
        return data['density'].in_cgs()

    def unitless_data(field, data):
            return np.ones(data['density'].shape)

    ds.add_field('density_alias_no_units', function=density_alias)
    ds.add_field('density_alias_unknown', function=density_alias,
                 units='unknown')
    ds.add_field('density_alias_wrong_units', function=density_alias,
                 units='m/s')

    assert_raises(YTFieldUnitError, get_data, ds, 'density_alias_no_units')
    assert_raises(YTFieldUnitError, get_data, ds, 'density_alias_wrong_units')

    dens = ad['density_alias_unknown']
    assert_equal(str(dens.units), 'g/cm**3')

    ds.add_field('dimensionless', function=unitless_data)
    ds.add_field('dimensionless_unknown', function=unitless_data,
                 units='unknown')
    ds.add_field('dimensionless_explicit', function=unitless_data, units='')
    ds.add_field('dimensionful', function=unitless_data, units='g/cm**3')

    assert_equal(str(ad['dimensionless'].units), 'dimensionless')
    assert_equal(str(ad['dimensionless_unknown'].units), 'dimensionless')
    assert_equal(str(ad['dimensionless_explicit'].units), 'dimensionless')
    assert_raises(YTFieldUnitError, get_data, ds, 'dimensionful')

if __name__ == "__main__":
    setup()
    for t in test_all_fields():
        t()
    test_add_deposited_particle_field()
    test_add_field_unit_semantics()
