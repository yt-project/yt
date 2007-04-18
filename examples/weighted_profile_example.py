import yt.lagos as lagos
import raven
import sys

a=lagos.EnzoHierarchy("DataDump0011")

v,c = a.findMax("Density")

region = lagos.EnzoSphere(a, c, 100./a["au"], ["CellMass","Density","Volume"])
bins,profiles=region.makeProfile(["Density", "NumberDensity"],10,1./a["au"],100./a["au"])

# Now we have the data, and we can do whatever with it.
