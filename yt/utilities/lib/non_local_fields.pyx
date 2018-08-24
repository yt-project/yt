"""
Cython tools for working with non local fields.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
cimport numpy as np

cimport cython

from libc.math cimport sqrt

from cykdtree.kdtree cimport PyKDTree, KDTree
from yt.utilities.lib.particle_kdtree_tools import knn_position
from yt.utilities.lib.bounded_priority_queue cimport BoundedPriorityQueue

def calculate_non_local_field(np.float64_t[:, ::1] output_buffer,
                              np.float64_t[:, :] particle_positions,
                              np.float64_t[:] mass,
                              np.float64_t[:] density,
                              np.float64_t[:] hsml,
                              np.float64_t[:, :] field_quantity,
                              PyKDTree kdtree,
                              spatial_type=None,
                              int num_neighbors=32,
                              kernel='cubic'):
    cdef int i, j, particle
    cdef np.float64_t prefactor
    cdef np.float64_t[::1] field_diff = np.zeros((1 + 2 * field_quantity.shape[1]))
    cdef np.float64_t[::1] kernel_values = np.zeros((3))
    cdef BoundedPriorityQueue queue = BoundedPriorityQueue(num_neighbors, True)

    # Get the function for the calculating the spatial derivatives
    spatial_type = divergence

    for i in range(particle_positions.shape[0]):
        # calculate the nearest neighbors here
        knn_position(particle_positions[i, :], particle_positions, queue,
                     kdtree)

        for j in range(num_neighbors):
            # Calculate the pre-factor and the difference in fields
            particle = queue.pids[j]
            prefactor = mass[particle] / denstiy[particle]
            field_diff[:] = field_quantity[particle, :] - field_quantity[i, :]

            # Calculate the kernel values here
            x = particle_positions[i, 0] - particle_positions[particle, 0]
            y = particle_positions[i, 1] - particle_positions[particle, 1]
            z = particle_positions[i, 2] - particle_positions[particle, 2]
            q = math.sqrt((x*x + y*y + z*z)) / hsml
            cubic_spline_grad(&kernel_values[0], q)

            # Add to the output buffer
            spatial_type(output_buffer[i, :], &kernel_values[0], &field_diff[0])

#TODO: fix this
cdef void cubic_spline_grad(np.float64_t * kernel_values,
                            np.float64_t x):
    cdef np.float64_t kernel
    # C is 8/pi
    cdef np.float64_t C = 2.5464790894703255
    if x <= 0.5:
        kernel = -12.0 * x + 18.0 * x * x
    elif x > 0.5 and x <= 1.0:
        kernel = -6.0 * (1.0 - x) * (1.0 - x)
    else:
        kernel = 0.0

    kernel_values[0] = kernel * C
    kernel_values[1] = kernel * C
    kernel_values[2] = kernel * C

# TODO: finish this
def output_vector(fields, spatial_type):
    '''
    This is a utility function to determine whether the output is a vector or a
    scalar quantity
    '''
    return True

cdef void gradient(np.float64_t * result, np.float64_t * kernel_values,
                   np.float64_t field):
    '''
    This takes in a scalar field and the 3 dimensional array of the gradient of
    the kernel values. This calculates the scalar field multiplied by the
    gradient of the kernel.
    '''
    result[0] += kernel_values[0] * field
    result[1] += kernel_values[1] * field
    result[2] += kernel_values[2] * field

cdef void divergence(np.float64_t * result, np.float64_t * kernel_values,
                     np.float64_t * field):
    '''
    This takes in a vector field and the array of the gradient of the kernel
    values. This calculates the dot product of the vector field and the gradient
    of the kernel.
    '''
    result[0] += kernel_values[0] * field[0]
    result[0] += kernel_values[1] * field[1]
    result[0] += kernel_values[2] * field[2]

cdef void curl(np.float64_t * result, np.float64_t * kernel_values,
               np.float64_t * field):
    '''
    This takes in a vector field and the array of the gradient of the kernel
    values. This calculates the cross product of the vector field and the
    gradient of the kernel.
    '''
    result[0] += field[1] * kernel_values[2] - field[2] * kernel_values[1]
    result[1] += field[2] * kernel_values[0] - field[0] * kernel_values[2]
    result[2] += field[0] * kernel_values[1] - field[1] * kernel_values[0]
