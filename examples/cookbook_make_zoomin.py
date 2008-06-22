from yt.mods import *

pf = get_pf() # Last argument on command line turned into an output file

n_frames = 200
min_dx = 250
frame_template = "frames/frame_%05i"

pc = raven.PlotCollection(pf)

for i in range(3):
    pl = pc.add_slice("Density",i)
    pl.add_callback(raven.UnitBoundaryCallback('pc'))

for i,v in enumerate(na.logspace(0,
             na.log10(pf.h.get_smallest_dx()*min_dx), n_frames)):
    pc.set_width(v,'1')
    fn=pc.save(frame_template % (i))
