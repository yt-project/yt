#
# This is a quick-and-dirty method of testing the parallel projection sketch
# I've created.  (Matt)
#
# It works very, very poorly.  There are numerous ways it could be improved.
# For starters, it does not even return a correct projection -- although making
# a correct one should only *add* to the time.  Unfortunately, the returns from
# this script are very poor; I'm not sure if this method is going to work.
# Perhaps a problem with the pickling and passing of arrays is to blame.  I
# guess we'll see.
#

from nws.sleigh import Sleigh
from nws.client import NetWorkSpace
import time
import numpy as na

from yt import ytcfg

ytcfg["lagos","serialize"] = "False"

NW = 16

s = Sleigh(workerCount = NW, nwsHost = "orange.slac.stanford.edu")
ws = NetWorkSpace('testing', serverHost = "orange.slac.stanford.edu")

s.eachWorker("import os\nttt=os.getcwd()\n")

initString = \
r"""
from nws.client import NetWorkSpace

ws = NetWorkSpace('testing')

import time

time.sleep(SleighRank)    # For scipy.weave
import yt.lagos as lagos
import yt.lagos.PointCombine as PointCombine

from yt import ytcfg

ytcfg["lagos","serialize"] = "False"
"""

s.eachWorker(initString)
#arr=s.eachWorker("v,c")

fn = "/nfs/slac/g/ki/ki02/mturk/work/07-09-04-ms3-restart/RedshiftOutput0032.dir/RedshiftOutput0032"

import yt.lagos as lagos
a = lagos.EnzoStaticOutput(fn)

from math import floor

time1a = time.time()
for level in range(a.h.maxLevel):
    #NG = a.h.numGrids
    grids = a.h.levelIndices[level]
    NG = len(grids)
    if NG < NW:
        for i in range(NW-1):
            ws.store('myGrids', [])
        ws.store('myGrids', grids)
    else:
        for i in range(NW):
            si = int(floor(float(NG)/NW * (i)))
            ei = int(floor(float(NG)/NW * (i+1)))
            #print i,na.mgrid[si:ei].min(), na.mgrid[si:ei].max()
            ws.store('myGrids', grids[na.mgrid[si:ei]])

    t1=time.time()
    s.eachWorker("lagos.DispatchedProjection('%s')" % fn)
    #lagos.DispatchedProjection(fn)
    t2=time.time()
    print "Level %s took %0.3e seconds" % (level, t2-t1)

    p=[]
    while 1:
        t = ws.fetchTry('myData')
        if not t: break
        p.append(t)
    #print p
time2a = time.time()

# Now we do it ourselves

time1b=time.time()
a.h.proj(0,"Density")
time2b=time.time()

print "Parallel : %0.3e" % (time2a-time1a)
print "Serial   : %0.3e" % (time2b-time1b)


s.stop()
