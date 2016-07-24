"""
Create smooth camera paths from keyframes.



"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import random
import numpy as np
from yt.visualization.volume_rendering.create_spline import create_spline

class Keyframes(object):
    def __init__(self, x, y, z=None, north_vectors=None, up_vectors=None,
                 times=None, niter=50000, init_temp=10.0, alpha=0.999,
                 fixed_start=False):
        r"""Keyframes for camera path generation.

        From a set of keyframes with position and optional up and
        north vectors, an interpolated camera path is generated.

        Parameters
        ----------
        x : array_like
            The x positions of the keyframes
        y : array_like
            The y positions of the keyframes
        z : array_like, optional
            The z positions of the keyframes. Default: 0.0
        north_vectors : array_like, optional
            The north vectors of the keyframes. Default: None
        up_vectors : array_like, optional
            The up vectors of the keyframes. Default: None
        times : array_like, optional
            The times of the keyframes. Default: arange(N)
        niter : integer, optional
            Maximum number of iterations to find solution. Default: 50000
        init_temp : float, optional
            Intital temperature for simulated annealing when finding a
            solution.  Lower initial temperatures result in an initial solution
            in first several iterations that changes more rapidly. Default: 10.0
        alpha : float, optional
            Exponent in cooling function in simulated annealing.  Must be < 1.
            In each iteration, the temperature_new = temperature_old * alpha.
            Default: 0.999
        fixed_start: boolean, optional
            If true, the first point never changes when searching for shortest
            path.  Default: False

        Examples
        --------

        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> from yt.visualization.volume_rendering.camera_path import *

        # Make a camera path from 10 random (x,y,z) keyframes
        >>> data = np.random.random.((10,3))
        >>> kf = Keyframes(data[:,0], data[:,1], data[:,2])
        >>> path = kf.create_path(250, shortest_path=False)

        # Plot the keyframes in the x-y plane and camera path
        plt.plot(kf.pos[:,0], kf.pos[:,1], 'ko')
        plt.plot(path['position'][:,0], path['position'][:,1])
        plt.savefig('path.png')
        """
        Nx = len(x)
        Ny = len(y)
        if z is not None:
            Nz = len(z)
            ndims = 3
        else:
            Nz = 1
            ndims = 2
        if Nx*Ny*Nz != Nx**ndims:
            print("Need Nx (%d) == Ny (%d) == Nz (%d)" % (Nx, Ny, Nz))
            raise RuntimeError
        self.nframes = Nx
        self.pos = np.zeros((Nx,3))
        self.pos[:,0] = x
        self.pos[:,1] = y
        if z is not None:
            self.pos[:,2] = z
        else:
            self.pos[:,2] = 0.0
        self.north_vectors = north_vectors
        self.up_vectors = up_vectors
        if times is None:
            self.times = np.arange(self.nframes)
        else:
            self.times = times
        self.cartesian_matrix()
        self.setup_tsp(niter, init_temp, alpha, fixed_start)

    def setup_tsp(self, niter=50000, init_temp=10.0, alpha=0.999,
                  fixed_start=False):
        r"""Setup parameters for Travelling Salesman Problem.

        Parameters
        ----------
        niter : integer, optional
            Maximum number of iterations to find solution. Default: 50000
        init_temp : float, optional
            Intital temperature for simulated annealing when finding a
            solution.  Lower initial temperatures result in an initial solution
            in first several iterations that changes more rapidly. Default: 10.0
        alpha : float, optional
            Exponent in cooling function in simulated annealing.  Must be < 1.
            In each iteration, the temperature_new = temperature_old * alpha.
            Default: 0.999
        fixed_start: boolean, optional
            If true, the first point never changes when searching for shortest
            path.  Default: False
        """
        # randomize tour
        self.tour = range(self.nframes)
        np.random.shuffle(self.tour)
        if fixed_start:
            first = self.tour.index(0)
            self.tour[0], self.tour[first] = self.tour[first], self.tour[0]
        self.max_iterations = niter
        self.initial_temp = init_temp
        self.alpha = alpha
        self.fixed_start = fixed_start
        self.best_score = None
        self.best = None

    def set_times(self, times):
        self.times = times
        
    def rand_seq(self):
        r"""
        Generates values in random order, equivlanet to using shuffle
        in random without generation all values at once.
        """
        values = range(self.nframes)
        for i in range(self.nframes):
            # pick a random index into remaining values
            j = i + int(random.random() * (self.nframes-i))
            # swap the values
            values[j], values[i] = values[i], values[j]
            # return the swapped value
            yield values[i]

    def all_pairs(self):
        r"""
        Generates all (i,j) pairs for (i,j) for 0-size
        """
        for i in self.rand_seq():
            for j in self.rand_seq():
                yield (i,j)
    
    def reversed_sections(self, tour):
        r"""
        Generator to return all possible variations where a section
        between two cities are swapped.
        """
        for i,j in self.all_pairs():
            if i == j: continue
            copy = tour[:]
            if i < j:
                copy[i:j+1] = reversed(tour[i:j+1])
            else:
                copy[i+1:] = reversed(tour[:j])
                copy[:j] = reversed(tour[i+1:])
            if self.fixed_start:
                ind = copy.index(0)
                copy[0], copy[ind] = copy[ind], copy[0]
            if copy != tour: # no point return the same tour
                yield copy

    def cartesian_matrix(self):
        r"""
        Create a distance matrix for the city coords that uses
        straight line distance
        """
        self.dist_matrix = np.zeros((self.nframes, self.nframes))
        xmat = np.zeros((self.nframes, self.nframes))
        xmat[:,:] = self.pos[:,0]
        dx = xmat - xmat.T
        ymat = np.zeros((self.nframes, self.nframes))
        ymat[:,:] = self.pos[:,1]
        dy = ymat - ymat.T
        zmat = np.zeros((self.nframes, self.nframes))
        zmat[:,:] = self.pos[:,2]
        dz = zmat - zmat.T
        self.dist_matrix = np.sqrt(dx*dx + dy*dy + dz*dz)

    def tour_length(self, tour):
        r"""
        Calculate the total length of the tour based on the distance
        matrix
        """
        total = 0
        num_cities = len(tour)
        for i in range(num_cities):
            j = (i+1) % num_cities
            city_i = tour[i]
            city_j = tour[j]
            total += self.dist_matrix[city_i, city_j]
        return -total

    def cooling(self):
        T = self.initial_temp
        while True:
            yield T
            T = self.alpha * T

    def prob(self, prev, next, temperature):
        if next > prev:
            return 1.0
        else:
            return np.exp( -abs(next-prev) / temperature )

    def get_shortest_path(self):
        r"""Determine shortest path between all keyframes.

        Parameters
        ----------
        None.
        """
        # this obviously doesn't work. When someone fixes it, remove the NOQA
        self.setup_tsp(niter, init_temp, alpha, fixed_start)  # NOQA
        num_eval = 1
        cooling_schedule = self.cooling()
        current = self.tour
        current_score = self.tour_length(current)
        for temperature in cooling_schedule:
            done = False
            # Examine moves around the current position
            for next in self.reversed_sections(current):
                if num_eval >= self.max_iterations:
                    done = True
                    break
                next_score = self.tour_length(next)
                num_eval += 1

                # Anneal.  Accept new solution if a random number is
                # greater than our "probability".
                p = self.prob(current_score, next_score, temperature)
                if random.random() < p:
                    current = next
                    self.current_score = next_score
                    if self.current_score > self.best_score:
                        #print num_eval, self.current_score, self.best_score, current
                        self.best_score = self.current_score
                        self.best = current
                    break

            if done:
                break
        self.pos = self.pos[self.tour,:]
        if self.north_vectors is not None:
            self.north_vectors = self.north_vectors[self.tour]
        if self.up_vectors is not None:
            self.up_vectors = self.up_vectors[self.tour]

    def create_path(self, npoints, path_time=None, tension=0.5, shortest_path=False):
        r"""Create a interpolated camera path from keyframes.

        Parameters
        ----------
        npoints : integer
            Number of points to interpolate from keyframes
        path_time : array_like, optional
            Times of interpolated points.  Default: Linearly spaced
        tension : float, optional
            Controls how sharp of a curve the spline takes.  A higher tension
            allows for more sharp turns.  Default: 0.5
        shortest_path : boolean, optional
            If true, estimate the shortest path between the keyframes.
            Default: False

        Returns
        -------
        path : dict
            Dictionary (time, position, north_vectors, up_vectors) of camera
            path.  Also saved to self.path.
        """
        self.npoints = npoints
        self.path = {"time": np.zeros(npoints),
                     "position": np.zeros((npoints, 3)),
                     "north_vectors": np.zeros((npoints,3)),
                     "up_vectors": np.zeros((npoints,3))}
        if shortest_path:
            self.get_shortest_path()
        if path_time is None:
            path_time = np.linspace(0, self.nframes, npoints)
        self.path["time"] = path_time
        for dim in range(3):
            self.path["position"][:,dim] = create_spline(self.times, self.pos[:,dim],
                                                         path_time, tension=tension)
            if self.north_vectors is not None:
                self.path["north_vectors"][:,dim] = \
                    create_spline(self.times, self.north_vectors[:,dim],
                                  path_time, tension=tension)
            if self.up_vectors is not None:
                self.path["up_vectors"][:,dim] = \
                    create_spline(self.times, self.up_vectors[:,dim],
                                  path_time, tension=tension)
        return self.path

    def write_path(self, filename="path.dat"):
        r"""Writes camera path to ASCII file

        Parameters
        ----------
        filename : string, optional
            Filename containing the camera path.  Default: path.dat
        """
        fp = open(filename, "w")
        fp.write("#%11s %12s %12s %12s %12s %12s %12s %12s %12s\n" %
                 ("x", "y", "z", "north_x", "north_y", "north_z",
                  "up_x", "up_y", "up_z"))
        for i in range(self.npoints):
            fp.write("%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n" %
                     (self.path["position"][i,0], self.path["position"][i,1],
                      self.path["position"][i,2], 
                      self.path["north_vectors"][i,0], self.path["north_vectors"][i,1],
                      self.path["north_vectors"][i,2], 
                      self.path["up_vectors"][i,0], self.path["up_vectors"][i,1],
                      self.path["up_vectors"][i,2]))
        fp.close()

