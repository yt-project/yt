import numpy as na
from yt.mods import *

def enhance(im, stdval=6.0, just_alpha=True):
    if just_alpha:
        nz = im[im>0.0]
        im[:] = im[:]/(nz.mean()+stdval*na.std(nz))
        im[im>1.0]=1.0
        im[im<0.0]=0.0
    else:
        for c in range(3):
            nz = im[:,:,c][im[:,:,c]>0.0]
            im[:,:,c] = im[:,:,c]/(nz.mean()+stdval*na.std(nz))
            del nz
        im[:,:][im>1.0]=1.0
        im[:,:][im<0.0]=0.0

if __name__ == 'main':
    im = na.zeros((256,256,3))
    line(im, 50,60,150,200)
    write_bitmap(im,'test_line.png')

