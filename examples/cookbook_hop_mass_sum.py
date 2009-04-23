from yt.mods import *

pf = get_pf() # last argument on the command line gets turned into an EnzoStaticOutput

hop_results = HaloFinder.HOPHaloFinder(pf, threshold=80.0)

def get_mass_results(hop_group):
    sphere = hop_group.get_sphere()
    baryon_mass = sphere["CellMassMsun"].sum()
    dm = sphere["creation_time"] < 0
    dm_mass = sphere["ParticleMassMsun"][dm].sum()
    star_mass = sphere["ParticleMassMsun"][~dm].sum()
    return "Total mass in HOP group %s is %0.5e (gas = %0.5e / dm = %0.5e / star = %0.5e)" % \
           (hop_group.id, baryon_mass + dm_mass + star_mass, baryon_mass, dm_mass, star_mass)

s = [get_mass_results(g) for g in hop_results]
print "\n".join(s)
