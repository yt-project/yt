"""
This file contains coordinate mappings between physical coordinates and those
defined on unit elements, as well as doing the corresponding intracell 
interpolation on finite element data.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from scipy.optimize import fsolve
cimport numpy as np
cimport cython

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


class P1Sampler2D:

    @staticmethod
    def map_real_to_unit(physical_coord, vertices):
    
        x = physical_coord[0]
        y = physical_coord[1]

        x1 = vertices[0, 0]
        y1 = vertices[0, 1]

        x2 = vertices[1, 0]
        y2 = vertices[1, 1]

        x3 = vertices[2, 0]
        y3 = vertices[2, 1]
    
        A = np.array([[1, x, y], [1, x1, y1], [1, x3, y3]])
        B = np.array([[1, x2, y2], [1, x1, y1], [1, x3, y3]])
        u = np.linalg.det(A) / np.linalg.det(B)

        C = np.array([[1, x, y], [1, x1, y1], [1, x2, y2]])
        D = np.array([[1, x3, y3], [1, x1, y1], [1, x2, y2]])
        v = np.linalg.det(C) / np.linalg.det(D)
    
        return np.array([u, v])

    @staticmethod
    def sample_at_unit_point(coord, vals):
        return vals[0]*(1 - coord[0] - coord[1]) + \
            vals[1]*coord[0] + vals[2]*coord[1]

    @classmethod
    def sample_at_real_point(cls, coord, vertices, vals):
        mapped_coord = cls.map_real_to_unit(coord, vertices)
        return cls.sample_at_unit_point(coord, vals)

class P1Sampler3D:

    @staticmethod
    def map_real_to_unit(physical_coord, vertices):
    
        x = physical_coord[0]
        y = physical_coord[1]
        z = physical_coord[2]

        x1 = vertices[0, 0]
        y1 = vertices[0, 1]
        z1 = vertices[0, 2]

        x2 = vertices[1, 0]
        y2 = vertices[1, 1]
        z2 = vertices[1, 2]
    
        x3 = vertices[2, 0]
        y3 = vertices[2, 1]
        z3 = vertices[2, 2]
    
        x4 = vertices[3, 0]
        y4 = vertices[3, 1]
        z4 = vertices[3, 2]
    
        b = np.array([x, y, z, 1])
        A = np.array([[x1, x2, x3, x4],
                      [y1, y2, y3, y4],
                      [z1, z2, z3, z4],
                      [1,  1,  1,  1] ])
    
        c = np.linalg.solve(A, b)
    
        return c

    @staticmethod
    def sample_at_unit_point(coord, vals):
        return vals[0]*coord[0] + vals[1]*coord[1] + \
            vals[2]*coord[2] + vals[3]*coord[3]

    @classmethod
    def sample_at_real_point(cls, coord, vertices, vals):
        mapped_coord = cls.map_real_to_unit(coord, vertices)
        return cls.sample_at_unit_point(coord, vals)


class Q1Sampler2D:

    def map_real_to_unit(self, physical_coord, vertices):
    
        # initial guess for the Newton solve
        x0 = np.array([0.0, 0.0])
        x = fsolve(self._f, x0, args=(vertices, physical_coord), fprime=self._J)
        return x

    @staticmethod
    def _f(x, v, phys_x):
        f1 = v[0][0]*(1-x[0])*(1-x[1]) + \
             v[1][0]*(1+x[0])*(1-x[1]) + \
             v[2][0]*(1-x[0])*(1+x[1]) + \
             v[3][0]*(1+x[0])*(1+x[1]) - 4.0*phys_x[0]
        f2 = v[0][1]*(1-x[0])*(1-x[1]) + \
             v[1][1]*(1+x[0])*(1-x[1]) + \
             v[2][1]*(1-x[0])*(1+x[1]) + \
             v[3][1]*(1+x[0])*(1+x[1]) - 4.0*phys_x[1]
        return np.array([f1, f2])

    @staticmethod
    def _J(x, v, phys_x):
        f11 = -(1-x[1])*v[0][0] + \
               (1-x[1])*v[1][0] - \
               (1+x[1])*v[2][0] + \
               (1+x[1])*v[3][0]
        f12 = -(1-x[0])*v[0][0] - \
               (1+x[0])*v[1][0] + \
               (1-x[0])*v[2][0] + \
               (1+x[0])*v[3][0]
        f21 = -(1-x[1])*v[0][1] + \
               (1-x[1])*v[1][1] - \
               (1+x[1])*v[2][1] + \
               (1+x[1])*v[3][1]
        f22 = -(1-x[0])*v[0][1] - \
               (1+x[0])*v[1][1] + \
               (1-x[0])*v[2][1] + \
               (1+x[0])*v[3][1]
        return np.array([[f11, f12], [f21, f22]])


class Q1Sampler3D:

    def map_real_to_unit(self, physical_coord, vertices):
    
        x0 = np.array([0.0, 0.0, 0.0])  # initial guess
        x = fsolve(self._f, x0, args=(vertices, physical_coord), fprime=self._J)
        return x

    @staticmethod
    def _f(x, v, phys_x):
        f0 = v[0][0]*(1-x[0])*(1-x[1])*(1-x[2]) + \
             v[1][0]*(1+x[0])*(1-x[1])*(1-x[2]) + \
             v[2][0]*(1-x[0])*(1+x[1])*(1-x[2]) + \
             v[3][0]*(1+x[0])*(1+x[1])*(1-x[2]) + \
             v[4][0]*(1-x[0])*(1-x[1])*(1+x[2]) + \
             v[5][0]*(1+x[0])*(1-x[1])*(1+x[2]) + \
             v[6][0]*(1-x[0])*(1+x[1])*(1+x[2]) + \
             v[7][0]*(1+x[0])*(1+x[1])*(1+x[2]) - 8.0*phys_x[0]
        f1 = v[0][1]*(1-x[0])*(1-x[1])*(1-x[2]) + \
             v[1][1]*(1+x[0])*(1-x[1])*(1-x[2]) + \
             v[2][1]*(1-x[0])*(1+x[1])*(1-x[2]) + \
             v[3][1]*(1+x[0])*(1+x[1])*(1-x[2]) + \
             v[4][1]*(1-x[0])*(1-x[1])*(1+x[2]) + \
             v[5][1]*(1+x[0])*(1-x[1])*(1+x[2]) + \
             v[6][1]*(1-x[0])*(1+x[1])*(1+x[2]) + \
             v[7][1]*(1+x[0])*(1+x[1])*(1+x[2]) - 8.0*phys_x[1]
        f2 = v[0][2]*(1-x[0])*(1-x[1])*(1-x[2]) + \
             v[1][2]*(1+x[0])*(1-x[1])*(1-x[2]) + \
             v[2][2]*(1-x[0])*(1+x[1])*(1-x[2]) + \
             v[3][2]*(1+x[0])*(1+x[1])*(1-x[2]) + \
             v[4][2]*(1-x[0])*(1-x[1])*(1+x[2]) + \
             v[5][2]*(1+x[0])*(1-x[1])*(1+x[2]) + \
             v[6][2]*(1-x[0])*(1+x[1])*(1+x[2]) + \
             v[7][2]*(1+x[0])*(1+x[1])*(1+x[2]) - 8.0*phys_x[2]
        return np.array([f0, f1, f2])

    @staticmethod
    def _J(x, v, phys_x):
    
        f00 = -(1-x[1])*(1-x[2])*v[0][0] + (1-x[1])*(1-x[2])*v[1][0] - \
               (1+x[1])*(1-x[2])*v[2][0] + (1+x[1])*(1-x[2])*v[3][0] - \
               (1-x[1])*(1+x[2])*v[4][0] + (1-x[1])*(1+x[2])*v[5][0] - \
               (1+x[1])*(1+x[2])*v[6][0] + (1+x[1])*(1+x[2])*v[7][0]
        f01 = -(1-x[0])*(1-x[2])*v[0][0] - (1+x[0])*(1-x[2])*v[1][0] + \
               (1-x[0])*(1-x[2])*v[2][0] + (1+x[0])*(1-x[2])*v[3][0] - \
               (1-x[0])*(1+x[2])*v[4][0] - (1+x[0])*(1+x[2])*v[5][0] + \
               (1-x[0])*(1+x[2])*v[6][0] + (1+x[0])*(1+x[2])*v[7][0]
        f02 = -(1-x[0])*(1-x[1])*v[0][0] - (1+x[0])*(1-x[1])*v[1][0] - \
               (1-x[0])*(1+x[1])*v[2][0] - (1+x[0])*(1+x[1])*v[3][0] + \
               (1-x[0])*(1-x[1])*v[4][0] + (1+x[0])*(1-x[1])*v[5][0] + \
               (1-x[0])*(1+x[1])*v[6][0] + (1+x[0])*(1+x[1])*v[7][0]
        

        f10 = -(1-x[1])*(1-x[2])*v[0][1] + (1-x[1])*(1-x[2])*v[1][1] - \
               (1+x[1])*(1-x[2])*v[2][1] + (1+x[1])*(1-x[2])*v[3][1] - \
               (1-x[1])*(1+x[2])*v[4][1] + (1-x[1])*(1+x[2])*v[5][1] - \
               (1+x[1])*(1+x[2])*v[6][1] + (1+x[1])*(1+x[2])*v[7][1]
        f11 = -(1-x[0])*(1-x[2])*v[0][1] - (1+x[0])*(1-x[2])*v[1][1] + \
               (1-x[0])*(1-x[2])*v[2][1] + (1+x[0])*(1-x[2])*v[3][1] - \
               (1-x[0])*(1+x[2])*v[4][1] - (1+x[0])*(1+x[2])*v[5][1] + \
               (1-x[0])*(1+x[2])*v[6][1] + (1+x[0])*(1+x[2])*v[7][1]
        f12 = -(1-x[0])*(1-x[1])*v[0][1] - (1+x[0])*(1-x[1])*v[1][1] - \
               (1-x[0])*(1+x[1])*v[2][1] - (1+x[0])*(1+x[1])*v[3][1] + \
               (1-x[0])*(1-x[1])*v[4][1] + (1+x[0])*(1-x[1])*v[5][1] + \
               (1-x[0])*(1+x[1])*v[6][1] + (1+x[0])*(1+x[1])*v[7][1]
        
        f20 = -(1-x[1])*(1-x[2])*v[0][2] + (1-x[1])*(1-x[2])*v[1][2] - \
               (1+x[1])*(1-x[2])*v[2][2] + (1+x[1])*(1-x[2])*v[3][2] - \
               (1-x[1])*(1+x[2])*v[4][2] + (1-x[1])*(1+x[2])*v[5][2] - \
               (1+x[1])*(1+x[2])*v[6][2] + (1+x[1])*(1+x[2])*v[7][2]
        f21 = -(1-x[0])*(1-x[2])*v[0][2] - (1+x[0])*(1-x[2])*v[1][2] + \
               (1-x[0])*(1-x[2])*v[2][2] + (1+x[0])*(1-x[2])*v[3][2] - \
               (1-x[0])*(1+x[2])*v[4][2] - (1+x[0])*(1+x[2])*v[5][2] + \
               (1-x[0])*(1+x[2])*v[6][2] + (1+x[0])*(1+x[2])*v[7][2]
        f22 = -(1-x[0])*(1-x[1])*v[0][2] - (1+x[0])*(1-x[1])*v[1][2] - \
               (1-x[0])*(1+x[1])*v[2][2] - (1+x[0])*(1+x[1])*v[3][2] + \
               (1-x[0])*(1-x[1])*v[4][2] + (1+x[0])*(1-x[1])*v[5][2] + \
               (1-x[0])*(1+x[1])*v[6][2] + (1+x[0])*(1+x[1])*v[7][2]

        return np.array([[f00, f01, f02], [f10, f11, f12], [f20, f21, f22]])
