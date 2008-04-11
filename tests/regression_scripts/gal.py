""" Original with Weave:
Doing contour 9 / 10 (4.56276e-14 9.91625e-14)
yt.lagos   INFO       2008-04-10 17:04:45,085 Getting field Density from 197
yt.lagos   INFO       2008-04-10 17:04:46,284 Contouring over 1061337 cells with 30115 candidates
yt.lagos   INFO       2008-04-10 17:05:17,880 Getting field tempContours from 197
yt.lagos   INFO       2008-04-10 17:05:18,305 Identified 5 contours between 4.56276e-14 and 2.15510e-13
yt.lagos   INFO       2008-04-10 17:05:18,312 Getting field Contours from 197
yt.lagos   INFO       2008-04-10 17:05:19,330 Getting field GridIndices from 197
        Contour id 0.0 has: 4.56297e-14 9.88512e-14  (2351) (4 grids, 11081.0 11087.0)
        Contour id 1.0 has: 4.58381e-14 5.41707e-14  (13) (1 grids, 11079.0 11079.0)
        Contour id 2.0 has: 4.60523e-14 4.79728e-14  (19) (2 grids, 11079.0 11082.0)
        Contour id 3.0 has: 4.56287e-14 2.15439e-13  (26063) (32 grids, 11040.0 11103.0)
        Contour id 4.0 has: 4.56339e-14 7.18032e-14  (1669) (3 grids, 11039.0 11051.0)
"""

""" New with C extension:
yt.lagos   INFO       2008-04-11 00:18:03,943 Getting field Density from 197
yt.lagos   INFO       2008-04-11 00:18:05,109 Contouring over 1061337 cells with 30115 candidates
yt.lagos   INFO       2008-04-11 00:18:36,243 Getting field tempContours from 197
yt.lagos   INFO       2008-04-11 00:18:36,626 Identified 5 contours between 4.56276e-14 and 2.15510e-13
yt.lagos   INFO       2008-04-11 00:18:36,633 Getting field Contours from 197
yt.lagos   INFO       2008-04-11 00:18:37,601 Getting field GridIndices from 197
        Contour id 0.0 has: 3.43716e-14 9.88512e-14  (2352) (4 grids, 11081.0 11087.0)
        Contour id 1.0 has: 4.58381e-14 5.41707e-14  (13) (1 grids, 11079.0 11079.0)
        Contour id 2.0 has: 4.60523e-14 4.79728e-14  (19) (2 grids, 11079.0 11082.0)
        Contour id 3.0 has: 4.56287e-14 2.15439e-13  (26063) (32 grids, 11040.0 11103.0)
        Contour id 4.0 has: 4.56339e-14 7.18032e-14  (1669) (3 grids, 11039.0 11051.0)
"""


import sys
sys.path.insert(0,"/Users/matthewturk/Development/yt/trunk")
from yt.config import ytcfg

ytcfg["yt","LogLevel"] = '20'
ytcfg["yt","logFile"] = "False"
ytcfg["yt","suppressStreamLogging"] = "True"

con_field = "Temperature"
con_field = "Density"

import yt.lagos as lagos
import yt.raven as raven
import numpy as na
import scipy.weave as weave
import scipy.weave.converters as converters

print lagos, raven

#fn = "/Users/matthewturk/Research/data/yt_test_data/galaxy1360.dir/galaxy1360"
fn = "/Users/matthewturk/Research/data/DataDir0036/DataDump0036"
a = lagos.EnzoStaticOutput(fn)

v,c = a.h.find_max("Density")

sp = a.h.sphere(c,20000.0/a['au'])
sp2 = a.h.sphere(c,200.0/a['au'])

#lagos.identify_contours(sp, "Density", 1e-27, 1e-24)
#lagos.identify_contours(sp, "Density", 1e-33, 1e-22)
nc = 10
cons = na.logspace(na.log10(sp2[con_field].min()*0.9),
                   na.log10(sp2[con_field].max()),nc+1)
do_plot = 1
k = []
k_bad = []
for i in [8]:#range(nc):
    [grid.clear_data() for grid in sp._grids]
    mi, ma = cons[i], cons[i+1]
    print "Doing contour %s / %s (%0.5e %0.5e)" % (i+1,nc,mi,ma)
    mq=lagos.identify_contours(sp, con_field, mi, sp[con_field].max())
    for cid in mq:
        sp["Contours"][mq[cid]] = cid
    k.append(na.unique(sp["Contours"]).size)
    my_con = sp["Contours"]
    sp._flush_data_to_grids("Contours",-1)
    for cid in na.unique(sp["Contours"])[1:]:
        cid_ind = na.where(sp["Contours"] == cid)
        print "\tContour id %s has: %0.5e %0.5e  (%s) (%s grids, %s %s)" % \
            (cid, sp[con_field][cid_ind].min(), sp[con_field][cid_ind].max(),
             sp["Density"][cid_ind].size,
             na.unique(sp["GridIndices"][cid_ind]).size,
             na.unique(sp["GridIndices"][cid_ind]).min(),
             na.unique(sp["GridIndices"][cid_ind]).max(),
             )
    #sp = a.h.sphere(c,20000.0/a['au'])
    sp.get_data("Contours",in_grids=True)
    if not na.all(sp["Contours"] == my_con):
        khdoihiueh0
    g_bad = 0
    for grid in sp._grids:
        grid['dx'] = na.ones(grid.ActiveDimensions)*grid.dx
        grid['dy'] = na.ones(grid.ActiveDimensions)*grid.dy
        grid['dz'] = na.ones(grid.ActiveDimensions)*grid.dz
        g_bad+=lagos.check_neighbors(grid)
    print "GRID:",g_bad
    k_bad.append(lagos.check_neighbors(sp))
    print "SOURCE",k_bad[-1]
    if do_plot:
        #c=na.array([sp["x"][ind][-1],sp["y"][ind][-1],sp["z"][ind][-1]])
        pc = raven.PlotCollection(a, center=c)
        pc.add_slice("Density",0)
        pc.add_slice("Density",1)
        pc.add_slice("Density",2)
        for field in ["Density","Contours"]:
            pc.switch_field(field)
            pc.set_width(10000,'au')
            pc.save("test_bds_%02i_10000au" % i)
            pc.set_width(2000,'au')
            pc.save("test_bds_%02i_02000au" % i)
            pc.set_width(1000,'au')
            pc.save("test_bds_%02i_01000au" % i)
            pc.set_width(100,'au')
            pc.save("test_bds_%02i_00100au" % i)
    del sp.data["Contours"]
print k
print k_bad
