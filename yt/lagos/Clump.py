import yt.lagos as lagos
import numpy as na

cl_parameters = {'minCells':64}

class Clump(object):
    children = None
    def __init__(self, data, parent, field):
        self.parent = parent
        self.data = data
        self.field = field
        self.min = self.data[field].min()
        self.max = self.data[field].max()
        self.isBound = None

    def find_children(self, min, max = None):
        if self.children is not None:
            print "Wiping out existing children clumps."
        self.children = []
        if max is None: max = self.max
        contour_info = lagos.identify_contours(self.data, self.field, min, max)
        for cid in contour_info:
            new_clump = self.data.extract_region(contour_info[cid])
            self.children.append(Clump(new_clump, self, self.field))

def find_clumps(clump, min, max, d_clump):
    print "Finding clumps: min: %e, max: %e, step: %f" % (min, max, d_clump)
    if min >= max: return
    clump.find_children(min)

    these_children = []
    for child in clump.children:
        if (len(child.data["CellMassMsun"]) >= cl_parameters['minCells']):
            these_children.append(child)
        else:
            print "Eliminating clump with %d cells." % len(child.data["CellMassMsun"])
    clump.children = these_children

    if (len(clump.children) == 1):
        find_clumps(clump, min*d_clump, max, d_clump)
    else:
        for child in clump.children:
            find_clumps(child, min*d_clump, max, d_clump)

def write_clump_hierarchy(clump,level,f_ptr):
    for q in range(level):
        f_ptr.write("\t")
    f_ptr.write("Clump at level %d:\n" % level)
    write_clump_info(clump,level,f_ptr)
    f_ptr.write("\n")
    f_ptr.flush()
    if ((clump.children is not None) and (len(clump.children) > 0)):
        for child in clump.children:
            write_clump_hierarchy(child,(level+1),f_ptr)

def write_clumps(clump,level,f_ptr):
    if ((clump.children is None) or (len(clump.children) == 0)):
        for q in range(level):
            f_ptr.write("\t")
        f_ptr.write("Clump:\n")
        write_clump_info(clump,level,f_ptr)
        f_ptr.write("\n")
        f_ptr.flush()
    if ((clump.children is not None) and (len(clump.children) > 0)):
        for child in clump.children:
            write_clumps(child,0,f_ptr)

def write_clump_finder_parameters(f_ptr):
    for par in cl_parameters.keys():
        f_ptr.write("# %s: %s\n" % (par,cl_parameters[par]))

def write_clump_info(clump,level,f_ptr):
    if clump.isBound is None:
        clump.isBound = clump.data.quantities["IsBound"]()
    for q in range(level):
        f_ptr.write("\t")
    f_ptr.write("Cells: %d\n" % len(clump.data["CellMassMsun"]))
    for q in range(level):
        f_ptr.write("\t")
    f_ptr.write("Mass: %.6e Msolar\n" % clump.data["CellMassMsun"].sum())
    for q in range(level):
        f_ptr.write("\t")
    f_ptr.write("Jeans Mass (vol-weighted): %.6e Msolar\n" % \
                    (clump.data.quantities["WeightedAverageQuantity"]("JeansMassMsun","CellVolume")))
    for q in range(level):
        f_ptr.write("\t")
    f_ptr.write("Jeans Mass (mass-weighted): %.6e Msolar\n" % \
                    (clump.data.quantities["WeightedAverageQuantity"]("JeansMassMsun","CellMassMsun")))
    for q in range(level):
        f_ptr.write("\t")
    f_ptr.write("Max grid level: %d\n" % clump.data["GridLevel"].max())
    for q in range(level):
        f_ptr.write("\t")
    f_ptr.write("Min number density: %.6e cm^-3\n" % clump.data["NumberDensity"].min())
    for q in range(level):
        f_ptr.write("\t")
    f_ptr.write("Max number density: %.6e cm^-3\n" % clump.data["NumberDensity"].max())
    for q in range(level):
        f_ptr.write("\t")
    f_ptr.write("Bound: %s\n" % clump.isBound)
