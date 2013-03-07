from yt.testing import *
import numpy as np
from yt.data_objects.field_info_container import \
    FieldInfo
import yt.data_objects.universal_fields
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"
    np.seterr(all = 'ignore')

_sample_parameters = dict(
    axis = 0,
    center = np.array((0.0, 0.0, 0.0)),
    bulk_velocity = np.array((0.0, 0.0, 0.0)),
    normal = np.array((0.0, 0.0, 1.0)),
    cp_x_vec = np.array((1.0, 0.0, 0.0)),
    cp_y_vec = np.array((0.0, 1.0, 0.0)),
    cp_z_vec = np.array((0.0, 0.0, 1.0)),
)

_base_fields = ["Density", "x-velocity", "y-velocity", "z-velocity"]

def realistic_pf(fields, nprocs):
    pf = fake_random_pf(16, fields = fields, nprocs = nprocs)
    pf.parameters["HydroMethod"] = "streaming"
    pf.parameters["EOSType"] = 1.0
    pf.parameters["EOSSoundSpeed"] = 1.0
    pf.conversion_factors["Time"] = 1.0
    pf.conversion_factors.update( dict((f, 1.0) for f in fields) )
    pf.gamma = 5.0/3.0
    pf.current_redshift = 0.0001
    pf.hubble_constant = 0.7
    for unit in mpc_conversion:
        pf.units[unit+'h'] = pf.units[unit]
        pf.units[unit+'cm'] = pf.units[unit]
        pf.units[unit+'hcm'] = pf.units[unit]
    return pf

class TestFieldAccess(object):
    description = None

    def __init__(self, field_name, nproc):
        # Note this should be a field name
        self.field_name = field_name
        self.description = "Accessing_%s_%s" % (field_name, nproc)
        self.nproc = nproc

    def __call__(self):
        field = FieldInfo[self.field_name]
        deps = field.get_dependencies()
        fields = deps.requested + _base_fields
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
        dd1.field_parameters.update(_sample_parameters)
        dd2.field_parameters.update(_sample_parameters)
        v1 = dd1[self.field_name]
        conv = field._convert_function(dd1) or 1.0
        if not needs_spatial:
            assert_equal(v1, conv*field._function(field, dd2))
        if not skip_grids:
            for g in pf.h.grids:
                g.field_parameters.update(_sample_parameters)
                conv = field._convert_function(g) or 1.0
                v1 = g[self.field_name]
                g.clear_data()
                g.field_parameters.update(_sample_parameters)
                assert_equal(v1, conv*field._function(field, g))

def test_all_fields():
    for field in FieldInfo:
        if field.startswith("CuttingPlane"): continue
        if field.startswith("particle"): continue
        if field.startswith("CIC"): continue
        if field.startswith("WeakLensingConvergence"): continue
        if FieldInfo[field].particle_type: continue
        for nproc in [1, 4, 8]:
            yield TestFieldAccess(field, nproc)
