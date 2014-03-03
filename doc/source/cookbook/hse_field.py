from yt.mods import *

# Define the components of the gravitational acceleration vector field by taking the
# gradient of the gravitational potential

def _Grav_Accel_x(field, data) :

    # We need to set up stencils

    sl_left = slice(None,-2,None)
    sl_right = slice(2,None,None)
    div_fac = 2.0

    dx = div_fac * data['dx'].flat[0]

    gx  = data["Grav_Potential"][sl_right,1:-1,1:-1]/dx
    gx -= data["Grav_Potential"][sl_left, 1:-1,1:-1]/dx

    new_field = np.zeros(data["Grav_Potential"].shape, dtype='float64')
    new_field[1:-1,1:-1,1:-1] = -gx

    return new_field

def _Grav_Accel_y(field, data) :

    # We need to set up stencils

    sl_left = slice(None,-2,None)
    sl_right = slice(2,None,None)
    div_fac = 2.0

    dy = div_fac * data['dy'].flat[0]

    gy  = data["Grav_Potential"][1:-1,sl_right,1:-1]/dy
    gy -= data["Grav_Potential"][1:-1,sl_left ,1:-1]/dy

    new_field = np.zeros(data["Grav_Potential"].shape, dtype='float64')
    new_field[1:-1,1:-1,1:-1] = -gy

    return new_field

def _Grav_Accel_z(field, data) :

    # We need to set up stencils

    sl_left = slice(None,-2,None)
    sl_right = slice(2,None,None)
    div_fac = 2.0

    dz = div_fac * data['dz'].flat[0]

    gz  = data["Grav_Potential"][1:-1,1:-1,sl_right]/dz
    gz -= data["Grav_Potential"][1:-1,1:-1,sl_left ]/dz

    new_field = np.zeros(data["Grav_Potential"].shape, dtype='float64')
    new_field[1:-1,1:-1,1:-1] = -gz

    return new_field

# Define the components of the pressure gradient field

def _Grad_Pressure_x(field, data) :

    # We need to set up stencils

    sl_left = slice(None,-2,None)
    sl_right = slice(2,None,None)
    div_fac = 2.0

    dx = div_fac * data['dx'].flat[0]

    px  = data["pressure"][sl_right,1:-1,1:-1]/dx
    px -= data["pressure"][sl_left, 1:-1,1:-1]/dx

    new_field = np.zeros(data["pressure"].shape, dtype='float64')
    new_field[1:-1,1:-1,1:-1] = px

    return new_field

def _Grad_Pressure_y(field, data) :

    # We need to set up stencils

    sl_left = slice(None,-2,None)
    sl_right = slice(2,None,None)
    div_fac = 2.0

    dy = div_fac * data['dy'].flat[0]

    py  = data["pressure"][1:-1,sl_right,1:-1]/dy
    py -= data["pressure"][1:-1,sl_left ,1:-1]/dy

    new_field = np.zeros(data["pressure"].shape, dtype='float64')
    new_field[1:-1,1:-1,1:-1] = py

    return new_field

def _Grad_Pressure_z(field, data) :

    # We need to set up stencils

    sl_left = slice(None,-2,None)
    sl_right = slice(2,None,None)
    div_fac = 2.0

    dz = div_fac * data['dz'].flat[0]

    pz  = data["pressure"][1:-1,1:-1,sl_right]/dz
    pz -= data["pressure"][1:-1,1:-1,sl_left ]/dz

    new_field = np.zeros(data["pressure"].shape, dtype='float64')
    new_field[1:-1,1:-1,1:-1] = pz

    return new_field

# Define the "degree of hydrostatic equilibrium" field

def _HSE(field, data) :

    gx = data["density"]*data["Grav_Accel_x"]
    gy = data["density"]*data["Grav_Accel_y"]
    gz = data["density"]*data["Grav_Accel_z"]

    hx = data["Grad_Pressure_x"] - gx
    hy = data["Grad_Pressure_y"] - gy
    hz = data["Grad_Pressure_z"] - gz

    h = np.sqrt((hx*hx+hy*hy+hz*hz)/(gx*gx+gy*gy+gz*gz))

    return h

# Now add the fields to the database

add_field("Grav_Accel_x", function=_Grav_Accel_x, take_log=False,
          validators=[ValidateSpatial(1,["Grav_Potential"])])

add_field("Grav_Accel_y", function=_Grav_Accel_y, take_log=False,
          validators=[ValidateSpatial(1,["Grav_Potential"])])

add_field("Grav_Accel_z", function=_Grav_Accel_z, take_log=False,
          validators=[ValidateSpatial(1,["Grav_Potential"])])

add_field("Grad_Pressure_x", function=_Grad_Pressure_x, take_log=False,
          validators=[ValidateSpatial(1,["pressure"])])

add_field("Grad_Pressure_y", function=_Grad_Pressure_y, take_log=False,
          validators=[ValidateSpatial(1,["pressure"])])

add_field("Grad_Pressure_z", function=_Grad_Pressure_z, take_log=False,
          validators=[ValidateSpatial(1,["pressure"])])

add_field("HSE", function=_HSE, take_log=False)

# Open two files, one at the beginning and the other at a later time when there's a
# lot of sloshing going on.

pfi = load("GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0000")
pff = load("GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0350")

# Sphere objects centered at the cluster potential minimum with a radius of 200 kpc

sphere_i = pfi.h.sphere(pfi.domain_center, (200, "kpc"))
sphere_f = pff.h.sphere(pff.domain_center, (200, "kpc"))

# Average "degree of hydrostatic equilibrium" in these spheres

hse_i = sphere_i.quantities["WeightedAverageQuantity"]("HSE", "cell_mass")
hse_f = sphere_f.quantities["WeightedAverageQuantity"]("HSE", "cell_mass")

print "Degree of hydrostatic equilibrium initially: ", hse_i
print "Degree of hydrostatic equilibrium later: ", hse_f

# Just for good measure, take slices through the center of the domain of the two files

slc_i = SlicePlot(pfi, 2, ["density","HSE"], center=pfi.domain_center, width=(1.0, "mpc"))
slc_f = SlicePlot(pff, 2, ["density","HSE"], center=pff.domain_center, width=(1.0, "mpc"))

slc_i.save("initial")
slc_f.save("final")
