from yt.mods import *

pf = load('RedshiftOutput0000')

h = yt.lagos.HaloFinding.chainHF(pf, dm_only=False)

