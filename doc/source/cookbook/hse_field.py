import numpy as np
import yt

from yt.fields.field_plugin_registry import \
    register_field_plugin
from yt.fields.fluid_fields import \
    setup_gradient_fields


# Define the components of the gravitational acceleration vector field by
# taking the gradient of the gravitational potential
@register_field_plugin
def setup_my_fields(registry, ftype="gas", slice_info=None):
    setup_gradient_fields(registry, (ftype, "gravitational_potential"),
                          "cm ** 2 / s ** 2", slice_info)

# Define the "degree of hydrostatic equilibrium" field


@yt.derived_field(name='HSE', units=None, take_log=False,
                  display_name='Hydrostatic Equilibrium')
def HSE(field, data):

    gx = data["density"] * data["gravitational_potential_gradient_x"]
    gy = data["density"] * data["gravitational_potential_gradient_y"]
    gz = data["density"] * data["gravitational_potential_gradient_z"]

    hx = data["pressure_gradient_x"] - gx
    hy = data["pressure_gradient_y"] - gy
    hz = data["pressure_gradient_z"] - gz

    h = np.sqrt((hx * hx + hy * hy + hz * hz) / (gx * gx + gy * gy + gz * gz))

    return h


# Open a dataset from when there's a lot of sloshing going on.

ds = yt.load("GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0350")

# gradient operator requires periodic boundaries.  This dataset has
# open boundary conditions.  We need to hack it for now (this will be fixed
# in future version of yt)
ds.periodicity = (True, True, True)

# Take a slice through the center of the domain
slc = yt.SlicePlot(ds, 2, ["density", "HSE"], width=(1, 'Mpc'))

slc.save("hse")
