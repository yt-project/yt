from yt.mods import *

pf = load("my_data") # load "my_data"

hop_results = HaloFinder(pf)

for hop_group in hop_results:
    sphere = hop_group.get_sphere()
    baryon_mass, particle_mass = sphere.quantities["TotalQuantity"](
            ["CellMassMsun", "ParticleMassMsun"], lazy_reader=True)
    print "Total mass in HOP group %s is %0.5e (gas = %0.5e / particles = %0.5e)" % \
            (hop_group.id, baryon_mass + particle_mass, baryon_mass, particle_mass)
