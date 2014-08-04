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
dimensionless = sympify(1)

base_dimensions = [mass, length, time, temperature, angle, dimensionless]

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
charge   = (energy * length)**Rational(1, 2)  # proper 1/2 power

electric_field = charge / length**2
magnetic_field = electric_field

solid_angle = angle * angle

derived_dimensions = [rate, velocity, acceleration, jerk, snap, crackle, pop, 
                      momentum, force, energy, power, charge, electric_field, 
                      magnetic_field, solid_angle, flux, specific_flux, volume,
                      area]

dimensions = base_dimensions + derived_dimensions
