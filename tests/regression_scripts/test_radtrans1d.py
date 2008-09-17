from yt.lagos import RTIntegrator as RT
import numpy as na
import pylab, time

size = 1024

i_s = 1.0
o_s = na.zeros(size, dtype='float')
e = na.ones(size, dtype='float') * 0
e[size/2:] = 0.0
a = na.ones(size, dtype='float')
#a[:,:,:size/2] = 0.0

dx = na.ones(size, dtype='float') * 1.0/size
t1 = time.time()
RT.Transfer1D(i_s, o_s, e, a, dx, 0, size)
t2 = time.time()
print "Took %0.3e seconds for a %s box" % (t2-t1, o_s.shape)

pylab.clf()
pylab.plot(range(size), o_s)
pylab.savefig("RT/line_plot_0.png")
