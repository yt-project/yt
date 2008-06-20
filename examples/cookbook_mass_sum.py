from yt.mods import *

pf = get_pf() # last argument on the command line gets turned into an EnzoStaticOutput

sp = pf.h.sphere([0.5,0.5,0.5], 1.0) # Everything, no pre-loading of fields
baryon_mass = sp["CellMassMsun"].sum()
dm_mass = sp["ParticleMassMsun"][sp["particle_type"] == 1].sum()
star_mass = sp["ParticleMassMsun"][sp["particle_type"] == 2].sum()

print "Total mass in grids in %s is %0.5e (gas = %0.5e / dm = %0.5e / star = %0.5e)" % \
            (pf, baryon_mass + dm_mass + star_mass, baryon_mass, dm_mass, star_mass)

