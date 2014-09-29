"""
Base dimensions


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from sympy import Symbol, sympify, Rational

mass = Symbol("(mass)", positive=True)
length = Symbol("(length)", positive=True)
time = Symbol("(time)", positive=True)
temperature = Symbol("(temperature)", positive=True)
angle = Symbol("(angle)", positive=True)
current_si = Symbol("(current_si)", positive=True)
dimensionless = sympify(1)

base_dimensions = [mass, length, time, temperature, angle, current_si,
                   dimensionless]

#
# Derived dimensions
#

rate = 1 / time

velocity     = length / time
acceleration = length / time**2
jerk         = length / time**3
snap         = length / time**4
crackle      = length / time**5
pop          = length / time**6

area     = length * length
volume   = area * length
momentum = mass * velocity
force    = mass * acceleration
energy   = force * length
power    = energy / time
flux     = power / area
specific_flux = flux / rate

# Gaussian electromagnetic units
charge   = (energy * length)**Rational(1, 2)  # proper 1/2 power
current  = charge / time
electric_field = charge / length**2
magnetic_field = electric_field

# SI electromagnetic units
charge_si = current_si * time
electric_field_si = force / charge_si
magnetic_field_si = electric_field_si / velocity

solid_angle = angle * angle

derived_dimensions = [rate, velocity, acceleration, jerk, snap, crackle, pop, 
                      momentum, force, energy, power, charge, electric_field, 
                      magnetic_field, solid_angle, flux, specific_flux, volume,
                      area, current, current_si, charge_si, electric_field_si,
                      magnetic_field_si]

dimensions = base_dimensions + derived_dimensions
