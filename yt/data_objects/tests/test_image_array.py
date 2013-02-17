from yt.testing import *
from yt.data_objects.image_array import ImageArray
import numpy as np

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"
    np.seterr(all = 'ignore')

def test_rgba_rescale():
    im = np.zeros([64,128,4])
    for i in xrange(im.shape[0]):
        for k in xrange(im.shape[2]):
            im[i,:,k] = np.linspace(0.,10.*k, im.shape[1])
    im_arr = ImageArray(im)

    new_im = im_arr.rescale(inline=False)
    yield assert_equal, im_arr[:,:,:3].max(), 2*10.
    yield assert_equal, im_arr[:,:,3].max(), 3*10.
    yield assert_equal, new_im[:,:,:3].sum(axis=2).max(), 1.0 
    yield assert_equal, new_im[:,:,3].max(), 1.0

    im_arr.rescale()
    yield assert_equal, im_arr[:,:,:3].sum(axis=2).max(), 1.0
    yield assert_equal, im_arr[:,:,3].max(), 1.0

def test_imarr_writepng():
    im = np.zeros([64,128,4])
    for i in xrange(im.shape[0]):
        for k in xrange(im.shape[2]):
            im[i,:,k] = np.linspace(0.,10.*k, im.shape[1])

    im_arr = ImageArray(im)
    im_arr.write_png('standard.png')
    im_arr.write_png('non-scaled.png', rescale=False)
    im_arr.write_png('black_bg.png', background='black')
    im_arr.write_png('white_bg.png', background='white')
    im_arr.write_png('green_bg.png', background=[0,1,0,1])
    im_arr.write_png('transparent_bg.png', background=None)

