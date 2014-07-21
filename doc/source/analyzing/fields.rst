Particle Fields
====================================
Naturally, particle fields contain properties of particles rather than
grid cells.  Many of these fields have corresponding grid fields that
can be populated by "depositing" the particle values onto a yt grid.

General Particle Fields
------------------------------------
Every particle will contain both a ``particle_position`` and ``particle_velocity``
that tracks the position and velocity (respectively) in code units.


SPH Fields
------------------------------------
For gas particles from SPH simulations, each particle will typically carry
a field for the smoothing length `h`, which is roughly equivalent to 
`(m/\rho)^{1/3}`, where `m` and `rho` are the particle mass and density 
respectively.  This can be useful for doing neighbour finding.

