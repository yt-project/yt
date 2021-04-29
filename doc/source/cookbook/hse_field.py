import numpy as np

import yt

# Open a dataset from when there's a lot of sloshing going on.

ds = yt.load("GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0350")

# Define the components of the gravitational acceleration vector field by
# taking the gradient of the gravitational potential
grad_fields = ds.add_gradient_fields(("gas", "gravitational_potential"))

# We don't need to do the same for the pressure field because yt already
# has pressure gradient fields. Now, define the "degree of hydrostatic
# equilibrium" field.


def _hse(field, data):
    # Remember that g is the negative of the potential gradient
    gx = -data[("gas", "density")] * data[("gas", "gravitational_potential_gradient_x")]
    gy = -data[("gas", "density")] * data[("gas", "gravitational_potential_gradient_y")]
    gz = -data[("gas", "density")] * data[("gas", "gravitational_potential_gradient_z")]
    hx = data[("gas", "pressure_gradient_x")] - gx
    hy = data[("gas", "pressure_gradient_y")] - gy
    hz = data[("gas", "pressure_gradient_z")] - gz
    h = np.sqrt((hx * hx + hy * hy + hz * hz) / (gx * gx + gy * gy + gz * gz))
    return h


ds.add_field(
    ("gas", "HSE"),
    function=_hse,
    units="",
    take_log=False,
    display_name="Hydrostatic Equilibrium",
    sampling_type="cell",
)

# The gradient operator requires periodic boundaries.  This dataset has
# open boundary conditions.
ds.force_periodicity()

# Take a slice through the center of the domain
slc = yt.SlicePlot(ds, 2, [("gas", "density"), ("gas", "HSE")], width=(1, "Mpc"))

slc.save("hse")
