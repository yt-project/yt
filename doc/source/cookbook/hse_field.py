### THIS RECIPE IS CURRENTLY BROKEN IN YT-3.0
### DO NOT TRUST THIS RECIPE UNTIL THIS LINE IS REMOVED

import numpy as np
import yt

# Define the components of the gravitational acceleration vector field by
# taking the gradient of the gravitational potential

@yt.derived_field(name='grav_accel_x', units='cm/s**2', take_log=False)
def grav_accel_x(field, data):

    # We need to set up stencils

    sl_left = slice(None, -2, None)
    sl_right = slice(2, None, None)
    div_fac = 2.0

    dx = div_fac * data['dx'].flat[0]

    gx = data["gravitational_potential"][sl_right, 1:-1, 1:-1]/dx
    gx -= data["gravitational_potential"][sl_left, 1:-1, 1:-1]/dx

    new_field = np.zeros(data["gravitational_potential"].shape,
                         dtype='float64')*gx.unit_array
    new_field[1:-1, 1:-1, 1:-1] = -gx

    return new_field


@yt.derived_field(name='grav_accel_y', units='cm/s**2', take_log=False)
def grav_accel_y(field, data):

    # We need to set up stencils

    sl_left = slice(None, -2, None)
    sl_right = slice(2, None, None)
    div_fac = 2.0

    dy = div_fac * data['dy'].flat[0]

    gy = data["gravitational_potential"][1:-1, sl_right, 1:-1]/dy
    gy -= data["gravitational_potential"][1:-1, sl_left, 1:-1]/dy

    new_field = np.zeros(data["gravitational_potential"].shape,
                         dtype='float64')*gx.unit_array
    new_field[1:-1, 1:-1, 1:-1] = -gy

    return new_field


@yt.derived_field(name='grav_accel_z', units='cm/s**2', take_log=False)
def grav_accel_z(field, data):

    # We need to set up stencils

    sl_left = slice(None, -2, None)
    sl_right = slice(2, None, None)
    div_fac = 2.0

    dz = div_fac * data['dz'].flat[0]

    gz = data["gravitational_potential"][1:-1, 1:-1, sl_right]/dz
    gz -= data["gravitational_potential"][1:-1, 1:-1, sl_left]/dz

    new_field = np.zeros(data["gravitational_potential"].shape,
                         dtype='float64')*gx.unit_array
    new_field[1:-1, 1:-1, 1:-1] = -gz

    return new_field


# Define the components of the pressure gradient field


@yt.derived_field(name='grad_pressure_x', units='g/(cm*s)**2', take_log=False)
def grad_pressure_x(field, data):

    # We need to set up stencils

    sl_left = slice(None, -2, None)
    sl_right = slice(2, None, None)
    div_fac = 2.0

    dx = div_fac * data['dx'].flat[0]

    px = data["pressure"][sl_right, 1:-1, 1:-1]/dx
    px -= data["pressure"][sl_left, 1:-1, 1:-1]/dx

    new_field = np.zeros(data["pressure"].shape, dtype='float64')*px.unit_array
    new_field[1:-1, 1:-1, 1:-1] = px

    return new_field


@yt.derived_field(name='grad_pressure_y', units='g/(cm*s)**2', take_log=False)
def grad_pressure_y(field, data):

    # We need to set up stencils

    sl_left = slice(None, -2, None)
    sl_right = slice(2, None, None)
    div_fac = 2.0

    dy = div_fac * data['dy'].flat[0]

    py = data["pressure"][1:-1, sl_right, 1:-1]/dy
    py -= data["pressure"][1:-1, sl_left, 1:-1]/dy

    new_field = np.zeros(data["pressure"].shape, dtype='float64')*px.unit_array
    new_field[1:-1, 1:-1, 1:-1] = py

    return new_field


@yt.derived_field(name='grad_pressure_z', units='g/(cm*s)**2', take_log=False)
def grad_pressure_z(field, data):

    # We need to set up stencils

    sl_left = slice(None, -2, None)
    sl_right = slice(2, None, None)
    div_fac = 2.0

    dz = div_fac * data['dz'].flat[0]

    pz = data["pressure"][1:-1, 1:-1, sl_right]/dz
    pz -= data["pressure"][1:-1, 1:-1, sl_left]/dz

    new_field = np.zeros(data["pressure"].shape, dtype='float64')*px.unit_array
    new_field[1:-1, 1:-1, 1:-1] = pz

    return new_field


# Define the "degree of hydrostatic equilibrium" field

@yt.derived_field(name='HSE', units=None, take_log=False)
def HSE(field, data):

    gx = data["density"]*data["Grav_Accel_x"]
    gy = data["density"]*data["Grav_Accel_y"]
    gz = data["density"]*data["Grav_Accel_z"]

    hx = data["Grad_Pressure_x"] - gx
    hy = data["Grad_Pressure_y"] - gy
    hz = data["Grad_Pressure_z"] - gz

    h = np.sqrt((hx*hx+hy*hy+hz*hz)/(gx*gx+gy*gy+gz*gz))*gx.unit_array

    return h


# Open two files, one at the beginning and the other at a later time when
# there's a lot of sloshing going on.

dsi = yt.load("GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0000")
dsf = yt.load("GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0350")

# Sphere objects centered at the cluster potential minimum with a radius
# of 200 kpc

sphere_i = dsi.sphere(dsi.domain_center, (200, "kpc"))
sphere_f = dsf.sphere(dsf.domain_center, (200, "kpc"))

# Average "degree of hydrostatic equilibrium" in these spheres

hse_i = sphere_i.quantities["WeightedAverageQuantity"]("HSE", "cell_mass")
hse_f = sphere_f.quantities["WeightedAverageQuantity"]("HSE", "cell_mass")

print "Degree of hydrostatic equilibrium initially: ", hse_i
print "Degree of hydrostatic equilibrium later: ", hse_f

# Just for good measure, take slices through the center of the domains
# of the two files

slc_i = yt.SlicePlot(dsi, 2, ["density", "HSE"], center=dsi.domain_center,
                     width=(1.0, "Mpc"))
slc_f = yt.SlicePlot(dsf, 2, ["density", "HSE"], center=dsf.domain_center,
                     width=(1.0, "Mpc"))

slc_i.save("initial")
slc_f.save("final")
