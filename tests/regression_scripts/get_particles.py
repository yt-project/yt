import sys
sys.path.insert(0,"/Users/matthewturk/Development/yt/trunk/")

import numpy as na
import yt.lagos as lagos

a = lagos.EnzoStaticOutput("galaxy1800.dir/galaxy1800")

sp = a.h.sphere([0.5,0.5,0.5], 1.0)

print sp["ParticleMassMsun"].max()
print sp["particle_velocity_x"].max()

print sp['particle_velocity_x'][sp['particle_type']==2].size
print sp['particle_velocity_x'][sp['particle_type']==1].size

for grid in a.h.grids:
    LE = grid.LeftEdge
    RE = grid.RightEdge
    center = (LE + RE)/2.0
    grid.set_field_parameter('center', center)
    rad = grid["RadiusCode"].max()
    reg = a.h.region(center, LE, RE)
    sph = a.h.sphere(center, rad)
    print "Checking region corresponding to grid", center, LE, RE
    # Now let's manually check our sizes
    pos = na.array([reg['particle_position_x'],
                    reg['particle_position_y'],
                    reg['particle_position_z']]).transpose()
    spos = na.array([sph['particle_position_x'],
                     sph['particle_position_y'],
                     sph['particle_position_z']]).transpose()
    mrad = na.array( na.sqrt(((spos - center)**2.0).sum(axis=1)))
    if na.any( na.all(pos > RE,axis=1) | na.all(pos < LE,axis=1) ) \
       or na.any(mrad > rad):
        print grid, na.any(mrad > rad), na.any( na.all(pos > RE,axis=1) | na.all(pos < LE,axis=1) )
        raise KeyError
