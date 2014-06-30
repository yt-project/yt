#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import numpy as np
from math import sin, cos

def identity():
	return np.matrix([[1.0, 0.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0, 0.0],
                      [0.0, 0.0, 1.0, 0.0],
                      [0.0, 0.0, 0.0, 1.0]]).astype(np.float32)
def standard():
	return np.matrix([[1.0, 0.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0, 0.0],
                      [0.0, 0.0, 1.0, 4.0],
                      [0.0, 0.0, -1.0, 0.0]]).astype(np.float32)

def perspective(eye):
    [ex, ey, ez] = eye
    a = -(ex/ez)
    b = -(ey/ez)
    c =  1.0/ez
    return np.matrix([[1.0, 0.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0, 0.0],
                      [  a,   b, 1.0,   c],
                      [0.0, 0.0, 0.0, 0.0]]).astype(np.float32)

def rotate(k, vec):
    [l, m, n] = vec
    return np.matrix([[l*l*(1.0-cos(k))+cos(k),
                       l*m*(1.0-cos(k))+n*sin(k),
                       l*n*(1.0-cos(k))-m*sin(k),
                       0.0],
                      [m*l*(1.0-cos(k))-n*sin(k),
                       m*m*(1.0-cos(k))+cos(k),
                       m*n*(1.0-cos(k))+l*sin(k),
                       0.0],
                      [n*l*(1.0-cos(k))+m*sin(k),
                       n*m*(1.0-cos(k))-l*sin(k),
                       n*n*(1.0-cos(k))+cos(k),
                       0.0],
                      [0.0, 0.0, 0.0, 1.0]]).astype(np.float32)

def rotate_x(k):
    return np.matrix([[1.0,    0.0,     0.0, 0.0],
                      [0.0, cos(k), -sin(k), 0.0],
                      [0.0, sin(k),  cos(k), 0.0],
                      [0.0,    0.0,     0.0, 1.0]]).astype(np.float32)

def rotate_y(k):
    return np.matrix([[ cos(k), 0.0, sin(k), 0.0],
                      [    0.0, 1.0,    0.0, 0.0],
                      [-sin(k), 0.0, cos(k), 0.0],
                      [    0.0, 0.0,    0.0, 1.0]]).astype(np.float32)

def rotate_z(k):
    return np.matrix([[cos(k), -sin(k), 0.0, 0.0],
                      [sin(k),  cos(k), 0.0, 0.0],
                      [   0.0,     0.0, 1.0, 0.0],
                      [   0.0,     0.0, 0.0, 1.0]]).astype(np.float32)

def scale(vec):
    [x, y, z] = vec
    return np.matrix([[x,   0.0, 0.0, 0.0],
                      [0.0,   y, 0.0, 0.0],
                      [0.0, 0.0,   z, 0.0],
                      [0.0, 0.0, 0.0, 1.0]]).astype(np.float32)

def translate(vec):
    [x, y, z] = vec
    return np.matrix([[1.0, 0.0, 0.0,   x],
                      [0.0, 1.0, 0.0,   y],
                      [0.0, 0.0, 1.0,   z],
                      [0.0, 0.0, 0.0, 1.0]]).astype(np.float32)
