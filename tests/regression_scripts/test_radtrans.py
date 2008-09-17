from yt.lagos import RTIntegrator as RT
import numpy as na
import pylab, time

size = 16

i_s = na.ones((size*2,size*2), dtype='float')
o_s = na.zeros((size,size,size), dtype='float')
e = na.ones((size,size,size), dtype='float') * 0
e[:,:,size/2:] = 0.0
a = na.ones((size,size,size), dtype='float')
#a[:,:,:size/2] = 0.0

def output_plots(image_field, axis):
    arr = image_field.swapaxes(0, axis)
    for i in range(arr.shape[axis]):
        print i
        pylab.clf()
        pylab.imshow(arr[i,...], interpolation='nearest')
        pylab.colorbar()
        pylab.clim(image_field.min(), image_field.max())
        pylab.savefig("RT/image_%s_%04i.png" % (axis,i))

dx = 1.0/size
t1 = time.time()
RT.Transfer3D(i_s, o_s, e, a, 2,4, 2,4, 2, 8, 2, 2, dx)
t2 = time.time()
print "Took %0.3e seconds for a %s box" % (t2-t1, o_s.shape)

pylab.clf()
pylab.plot(range(size), o_s[:,0,0].ravel())
pylab.savefig("RT/line_plot_0.png")

pylab.clf()
pylab.plot(range(size), o_s[0,:,0].ravel())
pylab.savefig("RT/line_plot_1.png")

pylab.clf()
pylab.plot(range(size), o_s[0,0,:].ravel())
pylab.savefig("RT/line_plot_2.png")

pylab.clf()
pylab.imshow(i_s, interpolation='nearest')
pylab.colorbar()
pylab.savefig("RT/resultant_image.png")

output_plots(o_s, 0)
output_plots(o_s, 1)
output_plots(o_s, 2)
