### THIS RECIPE IS CURRENTLY BROKEN IN YT-3.0
### DO NOT TRUST THIS RECIPE UNTIL THIS LINE IS REMOVED

import numpy as np
import yt

# Define the components of the gravitational acceleration vector field by
# taking the gradient of the gravitational potential

@yt.derived_field(name='gravitational_acceleration_x',
                  units='cm/s**2', take_log=False,
                  validators=[yt.ValidateSpatial(1,["gravitational_potential"])])
def gravitational_acceleration_x(field, data):

    # We need to set up stencils

    sl_left = slice(None, -2, None)
    sl_right = slice(2, None, None)
    div_fac = 2.0

    dx = div_fac * data['dx'][0]

    gx = data["gravitational_potential"][sl_right, 1:-1, 1:-1]/dx
    gx -= data["gravitational_potential"][sl_left, 1:-1, 1:-1]/dx

    new_field = np.zeros(data["gravitational_potential"].shape,
                         dtype='float64')*gx.uq
    new_field[1:-1, 1:-1, 1:-1] = -gx

    return new_field


@yt.derived_field(name='gravitational_acceleration_y',
                  units='cm/s**2', take_log=False,
                  validators=[yt.ValidateSpatial(1,["gravitational_potential"])])
def gravitational_acceleration_y(field, data):

    # We need to set up stencils

    sl_left = slice(None, -2, None)
    sl_right = slice(2, None, None)
    div_fac = 2.0

    dy = div_fac * data['dy'].flatten()[0]

    gy = data["gravitational_potential"][1:-1, sl_right, 1:-1]/dy
    gy -= data["gravitational_potential"][1:-1, sl_left, 1:-1]/dy

    new_field = np.zeros(data["gravitational_potential"].shape,
                         dtype='float64')*gy.uq

    new_field[1:-1, 1:-1, 1:-1] = -gy

    return new_field


@yt.derived_field(name='gravitational_acceleration_z',
                  units='cm/s**2', take_log=False,
                  validators=[yt.ValidateSpatial(1,["gravitational_potential"])])
def gravitational_acceleration_z(field, data):

    # We need to set up stencils

    sl_left = slice(None, -2, None)
    sl_right = slice(2, None, None)
    div_fac = 2.0

    dz = div_fac * data['dz'].flatten()[0]

    gz = data["gravitational_potential"][1:-1, 1:-1, sl_right]/dz
    gz -= data["gravitational_potential"][1:-1, 1:-1, sl_left]/dz

    new_field = np.zeros(data["gravitational_potential"].shape,
                         dtype='float64')*gz.uq
    new_field[1:-1, 1:-1, 1:-1] = -gz

    return new_field


# Define the components of the pressure gradient field


@yt.derived_field(name='grad_pressure_x', units='g/(cm*s)**2', take_log=False,
                  validators=[yt.ValidateSpatial(1,["pressure"])])
def grad_pressure_x(field, data):

    # We need to set up stencils

    sl_left = slice(None, -2, None)
    sl_right = slice(2, None, None)
    div_fac = 2.0

    dx = div_fac * data['dx'].flatten()[0]

    px = data["pressure"][sl_right, 1:-1, 1:-1]/dx
    px -= data["pressure"][sl_left, 1:-1, 1:-1]/dx

    new_field = np.zeros(data["pressure"].shape, dtype='float64')*px.uq
    new_field[1:-1, 1:-1, 1:-1] = px

    return new_field


@yt.derived_field(name='grad_pressure_y', units='g/(cm*s)**2', take_log=False,
                  validators=[yt.ValidateSpatial(1,["pressure"])])
def grad_pressure_y(field, data):

    # We need to set up stencils

    sl_left = slice(None, -2, None)
    sl_right = slice(2, None, None)
    div_fac = 2.0

    dy = div_fac * data['dy'].flatten()[0]

    py = data["pressure"][1:-1, sl_right, 1:-1]/dy
    py -= data["pressure"][1:-1, sl_left, 1:-1]/dy

    new_field = np.zeros(data["pressure"].shape, dtype='float64')*py.uq
    new_field[1:-1, 1:-1, 1:-1] = py

    return new_field


@yt.derived_field(name='grad_pressure_z', units='g/(cm*s)**2', take_log=False,
                  validators=[yt.ValidateSpatial(1,["pressure"])])
def grad_pressure_z(field, data):

    # We need to set up stencils

    sl_left = slice(None, -2, None)
    sl_right = slice(2, None, None)
    div_fac = 2.0

    dz = div_fac * data['dz'].flatten()[0]

    pz = data["pressure"][1:-1, 1:-1, sl_right]/dz
    pz -= data["pressure"][1:-1, 1:-1, sl_left]/dz

    new_field = np.zeros(data["pressure"].shape, dtype='float64')*pz.uq
    new_field[1:-1, 1:-1, 1:-1] = pz

    return new_field


# Define the "degree of hydrostatic equilibrium" field

@yt.derived_field(name='HSE', units=None, take_log=False,
                  display_name='Hydrostatic Equilibrium')
def HSE(field, data):

    gx = data["density"]*data["gravitational_acceleration_x"]
    gy = data["density"]*data["gravitational_acceleration_y"]
    gz = data["density"]*data["gravitational_acceleration_z"]

    hx = data["grad_pressure_x"] - gx
    hy = data["grad_pressure_y"] - gy
    hz = data["grad_pressure_z"] - gz

    h = np.sqrt((hx*hx+hy*hy+hz*hz)/(gx*gx+gy*gy+gz*gz))

    return h


# Open a dataset from when there's a lot of sloshing going on.

ds = yt.load("GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0350")


# Take a slice through the center of the domain
slc = yt.SlicePlot(ds, 2, ["density", "HSE"], width=(1, 'Mpc'))

slc.save("hse")
