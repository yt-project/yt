from yt.mods import *

pf = load("my_data") # load this data file

sp = pf.h.all_data()
baryon_mass = sp["CellMassMsun"].sum()
particle_mass = sphere["ParticleMassMsun"].sum()

print "Total mass in grids in %s is %0.5e (gas = %0.5e / particles = %0.5e)" % \
            (pf, baryon_mass + particle_mass, baryon_mass, particle_mass)

