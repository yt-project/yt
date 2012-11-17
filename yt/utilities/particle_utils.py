"""
Utilities for particles.

Author: John ZuHone <jzuhone@gmail.com>
Affiliation: NASA/Goddard Space Flight Center
Homepage: http://yt-project.org/
License:
Copyright (C) 2012 John ZuHone.  All Rights Reserved.

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

class MatchPointsToGrids(object) :

    def __init__(self, h, x, y, z) :

        self.num_points = len(x)
        
        self.x = x
        self.y = y
        self.z = z
        
        self.point_grids = -np.ones((self.num_points)).astype("int32")
        self.all_idxs = np.arange((self.num_points))

        self.root_grids = h.select_grids(0)

        self.counter = 0
        
    def __call__(self) :
        """
        Loops over every root grid and follows down each's subtree
        until it finds leaf grids, then assigns that grid id to the points
        in that grid
        """

        self.pbar = get_pbar("Finding Grids", self.num_points)
        
        for root_grid in self.root_grids :
            
            self.check_positions(root_grid)

        self.pbar.finish()
        
        return self.point_grids
                                                                                        
    def check_positions(self, grid) :
        """
        Recursive function to traverse up and down the hierarchy
        """
        
        if len(self.all_idxs) == 0 : return # We've run out of points

        if len(grid.Children) > 0 :

            # This grid has children

            for child in grid.Children :

                self.check_positions(child)

        else :

            # We hit a leaf grid

            local_idxs = grid.is_in_grid(self.x[self.all_idxs], self.y[self.all_idxs],
                                         self.z[self.all_idxs])
            
            if np.any(local_idxs) :

                self.point_grids[self.all_idxs[local_idxs]] = grid.id-grid._id_offset
                self.all_idxs = self.all_idxs[np.logical_not(local_idxs)]

                self.counter += local_idxs.sum()

                self.pbar.update(self.counter)
                
            return
                                                                                                                                   
