from yt.mods import *
pf = get_pf()

pc = raven.PlotCollection(pf)

fn = "%(bn)s_%(width)010i_%(unit)s" # template for image file names

widths = []
widths += [(i, "kpc") for i in [1]]
widths += [(i, "pc") for i in [1000, 100, 10, 1]]
widths += [(i, "au") for i in [1000, 100, 10, 1]]
widths += [(i, "rsun") for i in [1000, 100, 10, 1]]

def linked_save(pf, pc):
    for width, unit in widths:
        pc.set_width(width,unit)
        vmin = min([p.norm.vmin for p in pc.plots])
        vmax = max([p.norm.vmax for p in pc.plots])
        pc.set_zlim(vmin,vmax)
        d = {'bn':pf.basename, 'width':width, 'unit':unit}
        print pc.save(fn % d)

pc.add_slice("MachNumber", 0)
pc.add_slice("MachNumber", 1)
pc.add_slice("MachNumber", 2)
linked_save(pf, pc)
