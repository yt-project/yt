from yt.testing import *
import numpy as np
from yt.utilities.cosmology import \
     Cosmology
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion
from yt.frontends.stream.fields import \
    StreamFieldInfo
from yt.units.yt_array import YTArray

def setup():
    global base_pf
    # Make this super teeny tiny
    fields, units = [], []

    for fname, (code_units, aliases, dn) in StreamFieldInfo.known_other_fields:
        fields.append(("gas", fname))
        units.append(code_units)
    base_pf = fake_random_pf(4, fields = fields, units = units)
    base_pf.h
    base_pf.cosmological_simulation = 1
    base_pf.cosmology = Cosmology()
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"
    np.seterr(all = 'ignore')

def get_params(pf):
    return dict(
        axis = 0,
        center = YTArray((0.0, 0.0, 0.0), "cm",
            registry = pf.unit_registry),
        bulk_velocity = YTArray((0.0, 0.0, 0.0),
            "cm/s", registry = pf.unit_registry),
        normal = YTArray((0.0, 0.0, 1.0),
            "", registry = pf.unit_registry),
        cp_x_vec = YTArray((1.0, 0.0, 0.0),
            "", registry = pf.unit_registry),
        cp_y_vec = YTArray((0.0, 1.0, 0.0),
            "", registry = pf.unit_registry),
        cp_z_vec = YTArray((0.0, 0.0, 1.0),
            "", registry = pf.unit_registry),
        omega_baryon = 0.04,
        observer_redshift = 0.0,
        source_redshift = 3.0,
    )

_base_fields = (("gas", "density"),
                ("gas", "velocity_x"),
                ("gas", "velocity_y"),
                ("gas", "velocity_z"))

def realistic_pf(fields, nprocs):
    np.random.seed(int(0x4d3d3d3))
    units = [base_pf._get_field_info(*f).units for f in fields]
    fields = [_strip_ftype(f) for f in fields]
    pf = fake_random_pf(16, fields = fields, units = units,
                        nprocs = nprocs)
    pf.parameters["HydroMethod"] = "streaming"
    pf.parameters["EOSType"] = 1.0
    pf.parameters["EOSSoundSpeed"] = 1.0
    pf.conversion_factors["Time"] = 1.0
    pf.conversion_factors.update( dict((f, 1.0) for f in fields) )
    pf.gamma = 5.0/3.0
    pf.current_redshift = 0.0001
    pf.cosmological_simulation = 1
    pf.hubble_constant = 0.7
    pf.omega_matter = 0.27
    pf.omega_lambda = 0.73
    pf.cosmology = Cosmology(hubble_constant=pf.hubble_constant,
                             omega_matter=pf.omega_matter,
                             omega_lambda=pf.omega_lambda,
                             unit_registry=pf.unit_registry)
    return pf

def _strip_ftype(field):
    if not isinstance(field, tuple):
        return field
    elif field[0] == "all":
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
        if self.field_name in base_pf.field_list:
            # Don't know how to test this.  We need some way of having fields
            # that are fallbacks be tested, but we don't have that now.
            return
        field = base_pf._get_field_info(*self.field_name)
        deps = field.get_dependencies(pf = base_pf)
        fields = deps.requested + list(_base_fields)
        skip_grids = False
        needs_spatial = False
        for v in field.validators:
            f = getattr(v, "fields", None)
            if f: fields += f
            if getattr(v, "ghost_zones", 0) > 0:
                skip_grids = True
            if hasattr(v, "ghost_zones"):
                needs_spatial = True
        pf = realistic_pf(fields, self.nproc)
        # This gives unequal sized grids as well as subgrids
        dd1 = pf.h.all_data()
        dd2 = pf.h.all_data()
        sp = get_params(pf)
        dd1.field_parameters.update(sp)
        dd2.field_parameters.update(sp)
        v1 = dd1[self.field_name]
        # No more conversion checking
        if not field.particle_type:
            assert_equal(v1, dd1[self.field_name])
        if not needs_spatial:
            with field.unit_registry(dd2):
                res = field._function(field, dd2)
                res = dd2.apply_units(res, field.units)
            assert_array_almost_equal_nulp(v1, res, 4)
        if not skip_grids:
            for g in pf.index.grids:
                g.field_parameters.update(sp)
                v1 = g[self.field_name]
                g.clear_data()
                g.field_parameters.update(sp)
                for ax in 'xyz':
                    assert_array_equal(g[ax].shape, v1.shape)
                r1 = field._function(field, g)
                assert_array_equal(r1.shape, v1.shape)
                with field.unit_registry(g):
                    res = field._function(field, g)
                    assert_array_equal(v1.shape, res.shape)
                    res = g.apply_units(res, field.units)
                assert_array_almost_equal_nulp(v1, res, 4)

def test_all_fields():
    for field in sorted(base_pf.field_info):
        if not isinstance(field, types.TupleType):
            field = ("unknown", field)
        finfo = base_pf._get_field_info(*field)
        if isinstance(field, types.TupleType):
            fname = field[0]
        else:
            fname = field
        if field[0] == "deposit": continue
        if field[1].find("beta_p") > -1: continue
        if finfo.particle_type: continue
        for nproc in [1, 4, 8]:
            test_all_fields.__name__ = "%s_%s" % (field, nproc)
            yield TestFieldAccess(field, nproc)

if __name__ == "__main__":
    setup()
    for t in test_all_fields():
        t()
