from yt.mods import *

pf = get_pf() # last argument on the command line gets turned into an EnzoStaticOutput

full_sphere = pf.h.sphere([0.5,0.5,0.5], 1.0) # Everything, no pre-loading of fields
hop_results = lagos.hop.HopList(full_sphere, 80.0) # threshold = 80

def print_mass_results(id, sphere):
    baryon_mass = sphere["CellMassMsun"].sum()
    dm_mass = sphere["ParticleMassMsun"][sphere["particle_type"] == 1].sum()
    star_mass = sphere["ParticleMassMsun"][sphere["particle_type"] == 2].sum()
    print "Total mass in grids in %s is %0.5e (gas = %0.5e / dm = %0.5e / star = %0.5e)" % \
                (id, baryon_mass + dm_mass + star_mass, baryon_mass, dm_mass, star_mass)

for i,group in enumerate(hop_results):
    halo_sphere = group.get_sphere()
    print_mass_results(group.id, halo_sphere)
