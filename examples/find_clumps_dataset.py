# This is a wrapper for the clump finder in Clump.py.
# Arguments:
#     prefix: name of file prefix for text output
#     data_object: object over which contouring is performed (region or sphere).
#     field: data field over which contours are made (example: "Density" or "AveragedDensity").
#     step: contouring stepsize.  The field minimum is multiplied by this value each round of 
#           the clump finding.
# Written by Britton Smith

import sys, time
from yt.mods import *
from math import *

def find_clumps_dataset(prefix,data_object,field,step):
    t1 = time.time()

    c_min = 10**floor(log10(data_object[field].min()))
    c_max = 10**floor(log10(data_object[field].max())+1)

    master_clump = Clump(data_object, None, field)
    find_clumps(master_clump, c_min, c_max, step)

    t2=time.time()
    print "Took %0.3e seconds" % (t2-t1)

    f = open('%s_clump_hierarchy.txt' % prefix,'w')
    write_clump_hierarchy(master_clump,0,f)
    f.close()

    f = open('%s_clumps.txt' % prefix,'w')
    write_clumps(master_clump,0,f)
    f.close()

    return master_clump
