from yt.config import ytcfg

ytcfg["yt","time_functions"] = "True"

from yt.mods import *

yt_counters("Full Time")

yt_counters("yt Hierarchy")
pf = load('data0005')

pf.h
yt_counters("yt Hierarchy")

h = yt.lagos.HaloFinding.parallelHF(pf, threshold=160.0, safety=2.5, \
dm_only=False,resize=True, fancy_padding=True, rearrange=True)

yt_counters("Writing Data")
h.write_out('dist-chain.out')
h.write_particle_lists_txt("chain")
h.write_particle_lists("chain")
yt_counters("Writing Data")

yt_counters("Full Time")
