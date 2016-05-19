"""
Utilities for line integral convolution annotation



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Code originally from Scipy Cookbook (http://wiki.scipy.org/Cookbook/LineIntegralConvolution),
# with bug fixed which leads to crash when non equal-size vector field in two
# dimensions is provided.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import numpy as np
cimport numpy as np
cimport cython

@cython.cdivision(True)
cdef void _advance_2d(double vx, double vy,
                     int* x, int* y,
                     double* fx, double* fy,
                     int w, int h):
    cdef double tx, ty
    if vx>=0:
        tx = (1-fx[0])/vx
    else:
        tx = -fx[0]/vx
    if vy>=0:
        ty = (1-fy[0])/vy
    else:
        ty = -fy[0]/vy
    if tx<ty:
        if vx>=0:
            x[0]+=1
            fx[0]=0
        else:
            x[0]-=1
            fx[0]=1
        fy[0]+=tx*vy
    else:
        if vy>=0:
            y[0]+=1
            fy[0]=0
        else:
            y[0]-=1
            fy[0]=1
        fx[0]+=ty*vx
    if x[0]>=w:
        x[0]=w-1
    if x[0]<0:
        x[0]=0
    if y[0]<0:
        y[0]=0
    if y[0]>=h:
        y[0]=h-1

def line_integral_convolution_2d(
        np.ndarray[double, ndim=3] vectors,
        np.ndarray[double, ndim=2] texture,
        np.ndarray[double, ndim=1] kernel):
    cdef int i,j,l,x,y
    cdef int h,w,kernellen
    cdef double fx, fy
    cdef np.ndarray[double, ndim=2] result

    w = vectors.shape[0]
    h = vectors.shape[1]
    t = vectors.shape[2]

    kernellen = kernel.shape[0]
    result = np.zeros((w,h),dtype=np.double)

    vectors = vectors[...,::-1].copy()

    for i in range(w):
        for j in range(h):
            x = i
            y = j
            fx = 0.5
            fy = 0.5
            
            l = kernellen//2
            result[i,j] += kernel[l]*texture[x,y]
            while l<kernellen-1:
                _advance_2d(vectors[x,y,0],vectors[x,y,1],
                        &x, &y, &fx, &fy, w, h)
                l+=1
                result[i,j] += kernel[l]*texture[x,y]
            
            x = i
            y = j
            fx = 0.5
            fy = 0.5
            
            while l>0:
                _advance_2d(-vectors[x,y,0],-vectors[x,y,1],
                        &x, &y, &fx, &fy, w, h)
                l-=1
                result[i,j] += kernel[l]*texture[x,y]
                    
    return result
