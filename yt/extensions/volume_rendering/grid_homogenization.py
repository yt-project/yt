"""
A set of classes for homogenizing AMR data into fully-tiling, single-resolution
bricks.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/

License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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


import numpy as na

class GridFace(object):
    def __init__(self, grid, direction, left):
        self.direction = direction
        if left:
            self.coord = grid.LeftEdge[direction]
        else:
            self.coord = grid.RightEdge[direction]
        self.left_edge = grid.LeftEdge.copy()
        self.right_edge = grid.RightEdge.copy()
        self.left_edge[direction] = self.right_edge[direction] = self.coord

    def __and__(self, other):
        return self.intersect(other)

    def intersect(self, *args):
        """
        Returns a True or False based on intersection with either a set of
        (direction, coord, left_edge, right_edge) or another GridFace object.
        """
        if len(args) == 1:
            face = args[0]
            return self._intersect_coords(
                    face.direction, face.coord,
                    face.left_edge, face.right_edge)
        return self._intersect_coords(*args)

    def _intersect_coords(self, left_edge, right_edge):
        return (na.all(left_edge < self.right_edge) and
                na.all(right_edge > self.left_edge))

class GridFaces(object):
    def __init__(self, grids):
        self.faces = [ [], [], [] ]
        for grid in grids:
            for direction in range(3):
                self.faces[direction].append( GridFace(grid, direction, True) )
                self.faces[direction].append( GridFace(grid, direction, False) )
        for f in self.faces:
            f.sort(key = lambda a: a.coord)

    def __getitem__(self, item):
        return self.faces[item]

class ProtoPrism(object):
    def __init__(self, left_edge, right_edge, subgrid_faces):
        self.left_edge = left_edge
        self.right_edge = right_edge
        self.subgrid_faces = subgrid_faces
        # Not sure if we need these, or if we'll just tack on an extra subgrid
        # face for the domain boundaries.
        
    def sweep(self, direction = 0, stack = 0):
        # This is the sweep algorithm.  We sweep over a given direction, and if
        # we find an intersection, we have to sweep over another one in the
        # next direction.
        split_location = self.left_edge.copy()
        for face in self.subgrid_faces[direction]:
            split_location[direction] = face.coord
            if split_location[direction] <= self.left_edge[direction]:
                continue
            if split_location[direction] >= self.right_edge[direction]:
                return [self]
                # I thought this next line was necessary, but I am no longer
                # certain.  It would only be used if stack < 2.
                #return self.sweep((direction + 1) % 3, stack + 1)
            if face.intersect(self.left_edge, split_location):
                left, right = self.split(split_location)
                LC = left.sweep((direction + 1) % 3)
                RC = right.sweep(direction)
                return LC + RC
        raise RuntimeError

    def split(self, split_location):
        left = ProtoPrism(self.left_edge, split_location, self.subgrid_faces)
        right = ProtoPrism(split_location, self.right_edge, self.subgrid_faces)
        return (left, right)
