from yt.testing import *
from yt.data_objects.image_array import ImageArray
import numpy as np
import os
import tempfile
import shutil

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

def test_image_array_hdf5():
    # Perform I/O in safe place instead of yt main dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    im = np.zeros([64,128,3])
    for i in xrange(im.shape[0]):
        for k in xrange(im.shape[2]):
            im[i,:,k] = np.linspace(0.,0.3*k, im.shape[1])

    myinfo = {'field':'dinosaurs', 'east_vector':np.array([1.,0.,0.]), 
        'north_vector':np.array([0.,0.,1.]), 'normal_vector':np.array([0.,1.,0.]),  
        'width':0.245, 'units':'cm', 'type':'rendering'}

    im_arr = ImageArray(im, info=myinfo)
    im_arr.save('test_3d_ImageArray')

    im = np.zeros([64,128])
    for i in xrange(im.shape[0]):
        im[i,:] = np.linspace(0.,0.3*k, im.shape[1])

    myinfo = {'field':'dinosaurs', 'east_vector':np.array([1.,0.,0.]), 
        'north_vector':np.array([0.,0.,1.]), 'normal_vector':np.array([0.,1.,0.]),  
        'width':0.245, 'units':'cm', 'type':'rendering'}

    im_arr = ImageArray(im, info=myinfo)
    im_arr.save('test_2d_ImageArray')

    os.chdir(curdir)
    # clean up
    shutil.rmtree(tmpdir)

def test_image_array_rgb_png():
    # Perform I/O in safe place instead of yt main dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    im = np.zeros([64,128,3])
    for i in xrange(im.shape[0]):
        for k in xrange(im.shape[2]):
            im[i,:,k] = np.linspace(0.,10.*k, im.shape[1])

    im_arr = ImageArray(im)
    im_arr.write_png('standard.png')

def test_image_array_rgba_png():
    # Perform I/O in safe place instead of yt main dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    im = np.zeros([64,128,4])
    for i in xrange(im.shape[0]):
        for k in xrange(im.shape[2]):
            im[i,:,k] = np.linspace(0.,10.*k, im.shape[1])

    im_arr = ImageArray(im)
    im_arr.write_png('standard.png')
    im_arr.write_png('non-scaled.png', rescale=False)
    im_arr.write_png('black_bg.png', background='black')
    im_arr.write_png('white_bg.png', background='white')
    im_arr.write_png('green_bg.png', background=[0.,1.,0.,1.])
    im_arr.write_png('transparent_bg.png', background=None)


def test_image_array_background():
    # Perform I/O in safe place instead of yt main dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    im = np.zeros([64,128,4])
    for i in xrange(im.shape[0]):
        for k in xrange(im.shape[2]):
            im[i,:,k] = np.linspace(0.,10.*k, im.shape[1])

    im_arr = ImageArray(im)
    im_arr.rescale()
    new_im = im_arr.add_background_color([1.,0.,0.,1.], inline=False)
    new_im.write_png('red_bg.png')
    im_arr.add_background_color('black')
    im_arr.write_png('black_bg2.png')
 
    os.chdir(curdir)
    # clean up
    shutil.rmtree(tmpdir)













