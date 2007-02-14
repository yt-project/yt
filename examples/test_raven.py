from yt.lagos import *
import yt.raven as raven

h=EnzoHierarchy("/usr/work/mturk/07-01-22-fixed-gamma2/DataDump0006.dir/DataDump0006")

plot=raven.EnzoHippo(h)
plot.addSlice("Density",0)
