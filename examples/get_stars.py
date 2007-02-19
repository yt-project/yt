import yt.lagos

a = yt.lagos.EnzoHierarchy("galaxy0003")
a.exportParticlesPB("dm.pb", scale=100.0, filter=1)
a.exportParticlesPB("stars.pb", filter=2, scale=100.0, fields=["metallicity_fraction", "particle mass"])
