"""
AMR kD-Tree Tools

Authors: Samuel Skillman <samskillman@gmail.com>
Affiliation: University of Colorado at Boulder

Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Samuel Skillman.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import numpy as np
from yt.funcs import *


def receive_and_reduce(comm, incoming_rank, image, add_to_front):
    mylog.debug('Receiving image from %04i' % incoming_rank)
    #mylog.debug( '%04i receiving image from %04i'%(self.comm.rank,back.owner))
    arr2 = comm.recv_array(incoming_rank, incoming_rank).reshape(
        (image.shape[0], image.shape[1], image.shape[2]))

    if add_to_front:
        front = arr2
        back = image
    else:
        front = image
        back = arr2

    if image.shape[2] == 3:
        # Assume Projection Camera, Add
        np.add(image, front, image)
        return image

    ta = 1.0 - front[:, :, 3]
    np.maximum(ta, 0.0, ta)
    # This now does the following calculation, but in a memory
    # conservative fashion
    # image[:,:,i  ] = front[:,:,i] + ta*back[:,:,i]
    image = back.copy()
    for i in range(4):
        np.multiply(image[:, :, i], ta, image[:, :, i])
    np.add(image, front, image)
    return image


def send_to_parent(comm, outgoing_rank, image):
    mylog.debug('Sending image to %04i' % outgoing_rank)
    comm.send_array(image, outgoing_rank, tag=comm.rank)


def scatter_image(comm, root, image):
    mylog.debug('Scattering from %04i' % root)
    image = comm.mpi_bcast(image, root=root)
    return image
