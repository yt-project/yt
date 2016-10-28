import numpy as np

from yt import \
    load
from yt.testing import \
    fake_amr_ds, \
    fake_random_ds, \
    fake_particle_ds, \
    assert_almost_equal, \
    assert_equal, \
    assert_array_almost_equal_nulp, \
    assert_array_equal, \
    assert_raises, \
    requires_file
from yt.utilities.cosmology import \
    Cosmology
from yt.frontends.stream.fields import \
    StreamFieldInfo
from yt.units.yt_array import \
    array_like_field, \
    YTArray, YTQuantity
from yt.utilities.exceptions import \
    YTFieldUnitError, \
    YTFieldUnitParseError, \
    YTDimensionalityError

base_ds = None


def setup():
    global base_ds
    # Make this super teeny tiny
    fields, units = [], []

    for fname, (code_units, aliases, dn) in StreamFieldInfo.known_other_fields:
        fields.append(("gas", fname))
        units.append(code_units)

    base_ds = fake_random_ds(4, fields=fields, units=units, particles=20)

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
    global base_ds
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


class TestFieldAccess(object):
    description = None

    def __init__(self, field_name, nproc):
        # Note this should be a field name
        self.field_name = field_name
        self.description = "Accessing_%s_%s" % (field_name, nproc)
        self.nproc = nproc

    def __call__(self):
        global base_ds
        field = base_ds._get_field_info(*self.field_name)
        deps = field.get_dependencies(ds = base_ds)
        requested = deps.requested
        particle_fields = \
            ['particle_position_x', 'particle_position_y', 'particle_position_z',
             'particle_velocity_x', 'particle_velocity_y', 'particle_velocity_z',
             'particle_mass']
        fields = list(_base_fields)

        for rf in requested:
            if rf[0] == 'io' or rf[0] == 'all':
                if rf not in particle_fields or rf[1] not in particle_fields:
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
    global base_ds
    for field in sorted(base_ds.field_info):
        if field[1].find("beta_p") > -1:
            continue
        if field[1].find("vertex") > -1:
            # don't test the vertex fields for now
            continue
        if field in base_ds.field_list:
            # Don't know how to test this.  We need some way of having fields
            # that are fallbacks be tested, but we don't have that now.
            continue

        for nproc in [1, 4, 8]:
            test_all_fields.__name__ = "%s_%s" % (field, nproc)
            yield TestFieldAccess(field, nproc)

def test_add_deposited_particle_field():
    # NOT tested: "std", "mesh_id", "nearest" and "simple_smooth"
    global base_ds
    ad = base_ds.all_data()

    # Test "count", "sum" and "cic" method
    for method in ["count", "sum", "cic"]:
        fn = base_ds.add_deposited_particle_field(('io', 'particle_mass'), method)
        expected_fn = 'io_%s' if method == "count" else 'io_%s_mass'
        assert_equal(fn, ('deposit', expected_fn % method))
        ret = ad[fn]
        if method == "count":
            assert_equal(ret.sum(), ad['particle_ones'].sum())
        else:
            assert_almost_equal(ret.sum(), ad['particle_mass'].sum())

    # Test "weighted_mean" method
    fn = base_ds.add_deposited_particle_field(('io', 'particle_ones'), 'weighted_mean',
                                              weight_field='particle_ones')
    assert_equal(fn, ('deposit', 'io_avg_ones'))
    ret = ad[fn]
    # The sum should equal the number of cells that have particles
    assert_equal(ret.sum(), np.count_nonzero(ad[("deposit", "io_count")]))

@requires_file('GadgetDiskGalaxy/snapshot_200.hdf5')
def test_add_smoothed_particle_field():
    ds = load('GadgetDiskGalaxy/snapshot_200.hdf5')
    fn = ds.add_smoothed_particle_field(('PartType0', 'particle_ones'))
    assert_equal(fn, ('deposit', 'PartType0_smoothed_particle_ones'))
    dd = ds.sphere('center', (500, 'code_length'))
    ret = dd[fn]
    assert_almost_equal(ret.sum(), 638.5652315154682)

def test_add_gradient_fields():
    global base_ds
    gfields = base_ds.add_gradient_fields(("gas","density"))
    gfields += base_ds.add_gradient_fields(("index", "ones"))
    field_list = [('gas', 'density_gradient_x'),
                  ('gas', 'density_gradient_y'),
                  ('gas', 'density_gradient_z'),
                  ('gas', 'density_gradient_magnitude'),
                  ('index', 'ones_gradient_x'),
                  ('index', 'ones_gradient_y'),
                  ('index', 'ones_gradient_z'),
                  ('index', 'ones_gradient_magnitude')]
    assert_equal(gfields, field_list)
    ad = base_ds.all_data()
    for field in field_list:
        ret = ad[field]
        if field[0] == 'gas':
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
        return data['density'].in_cgs()

    def unitless_data(field, data):
            return np.ones(data['density'].shape)

    ds.add_field(('gas','density_alias_no_units'), function=density_alias)
    ds.add_field(('gas','density_alias_auto'), function=density_alias,
                 units='auto', dimensions='density')
    ds.add_field(('gas','density_alias_wrong_units'), function=density_alias,
                 units='m/s')
    ds.add_field(('gas','density_alias_unparseable_units'), function=density_alias,
                 units='dragons')
    ds.add_field(('gas','density_alias_auto_wrong_dims'), function=density_alias,
                 units='auto', dimensions="temperature")
    assert_raises(YTFieldUnitError, get_data, ds, 'density_alias_no_units')
    assert_raises(YTFieldUnitError, get_data, ds, 'density_alias_wrong_units')
    assert_raises(YTFieldUnitParseError, get_data, ds,
                  'density_alias_unparseable_units')
    assert_raises(YTDimensionalityError, get_data, ds, 'density_alias_auto_wrong_dims')

    dens = ad['density_alias_auto']
    assert_equal(str(dens.units), 'g/cm**3')

    ds.add_field(('gas','dimensionless'), function=unitless_data)
    ds.add_field(('gas','dimensionless_auto'), function=unitless_data,
                 units='auto', dimensions='dimensionless')
    ds.add_field(('gas','dimensionless_explicit'), function=unitless_data, units='')
    ds.add_field(('gas','dimensionful'), function=unitless_data, units='g/cm**3')

    assert_equal(str(ad['dimensionless'].units), 'dimensionless')
    assert_equal(str(ad['dimensionless_auto'].units), 'dimensionless')
    assert_equal(str(ad['dimensionless_explicit'].units), 'dimensionless')
    assert_raises(YTFieldUnitError, get_data, ds, 'dimensionful')

def test_array_like_field():
    ds = fake_random_ds(4, particles=64)
    ad = ds.all_data()
    u1 = ad["particle_mass"].units
    u2 = array_like_field(ad, 1., ("all", "particle_mass")).units
    assert u1 == u2

def test_add_field_string():
    ds = fake_random_ds(16)
    ad = ds.all_data()

    def density_alias(field, data):
        return data['density']

    ds.add_field('density_alias', function=density_alias, units='g/cm**3')

    ad['density_alias']
    assert ds.derived_field_list[0] == 'density_alias'

def test_add_field_string_aliasing():
    ds = fake_random_ds(16)

    def density_alias(field, data):
        return data['density']

    ds.add_field('density_alias', function=density_alias, units='g/cm**3')

    ds.field_info['density_alias']
    ds.field_info['gas', 'density_alias']

    ds = fake_particle_ds()
    
    def pmass_alias(field, data):
        return data['particle_mass']
        
    ds.add_field('particle_mass_alias', function=pmass_alias, 
                 units='g', particle_type=True)

    ds.field_info['particle_mass_alias']
    ds.field_info['all', 'particle_mass_alias']
    

def test_morton_index():
    ds = fake_amr_ds()
    mi = ds.r["index", "morton_index"]
    mi2 = mi.view("uint64")
    assert_equal(np.unique(mi2).size, mi2.size)
    a1 = np.argsort(mi)
    a2 = np.argsort(mi2)
    assert_array_equal(a1, a2)

if __name__ == "__main__":
    setup()
    for t in test_all_fields():
        t()
    test_add_deposited_particle_field()
    test_add_field_unit_semantics()
    test_array_like_field()
    test_add_field_string()
