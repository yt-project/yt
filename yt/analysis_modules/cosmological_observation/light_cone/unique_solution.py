"""
Functions to generate unique light cone solutions.

Author: Britton Smith <brittons@origins.colorado.edu>
Affiliation: CASA/University of CO, Boulder
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Britton Smith.  All Rights Reserved.

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

import copy
import numpy as np
import random as rand
import sys

from yt.funcs import *
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_root_only

from .light_cone import \
     LightCone
from .common_n_volume import \
     common_volume

def project_unique_light_cones(lightcone, seed_file, field, **kwargs):
    r"""Make light cone projections using list of random
    seeds in a file.

    Parameters
    ----------
    lightcone : LightCone
        A LightCone object.
    seed_file : str
        File containing random seeds for light cone projections.
    field : str
        The field to be projected.

    Any additional parameters to project_light_cone should be
    given here.
    """

    seedList = _read_seed_file(seed_file)

    prefix = lightcone.output_prefix
    lightcone.calculate_light_cone_solution(seed=0)
    lastSeed = None

    for seed in seedList:
        if (seed['master'] != lastSeed):
            lightcone.rerandomize_light_cone_solution(seed['master'],
                                                      recycle=False)
            lastSeed = seed['master']
        if (seed['recycle'] is not None):
            lightcone.rerandomize_light_cone_solution(seed['recycle'],
                                                      recycle=True)

        lightcone.output_prefix = "%s_%s_%s" % \
          (prefix, seed['master'], seed['recycle'])
        lightcone.project_light_cone(field, **kwargs)

def find_unique_solutions(lightcone1, solutions=10, seed=None,
                          max_overlap=0.25, failures=10,
                          recycle=True, filename='unique.dat'):
    r"""Find a set of random seeds that will give light cones
    will minimal volume overlap.

    Parameters
    ----------
    lightcone1 : LightCone
        A LightCone object.
    solutions : int
        The desired number of unique solutions.
        Default: 10
    seed : int
        The random seed for generating light cone randomizations.
        Default: None
    max_overlap : float
        The maximum fractional volume an accepted solution can have
        in common with the other solutions.
        Default: 0.25
    failures : int
        The number of failed attempts before ending.
        Default: 10
    recycle : bool
        If True, attempt to use existing realizations.  This
        will speed up light cone projectinos.
        Default: True
    filename : str
        Output file name for the list of unique seeds.
        Default : 'unique.dat'

    """

    lightcone2 = copy.deepcopy(lightcone1)
    lightcone1.calculate_light_cone_solution(seed=0)
    lightcone2.calculate_light_cone_solution(seed=0)

    unique_seeds = []
    if recycle:
        master = None
    new_recycle_seed = None
    fails = 0
    recycle_fails = 0

    max_common = 0.0

    # Need to continuall save and reset the state of the random
    # number generator since it is being reset by the light
    # cone generator.
    if seed is None:
        state = None
    else:
        rand.seed(seed)
        state = rand.getstate()

    fail_digits = str(int(np.log10(failures))+1)

    while (len(unique_seeds) < solutions):
        # Create new random seed.
        if (recycle and master is not None):
            new_seed = master
            if state is not None: rand.setstate(state)
            new_recycle_seed = rand.randint(1, 1e9)
            state = rand.getstate()
        else:
            if state is not None: rand.setstate(state)
            new_seed = rand.randint(1, 1e9)
            state = rand.getstate()
            if recycle:
                master = new_seed
                recycle_fails = 0
            new_recycle_seed = None

        sys.stderr.write(("Unique solutions: %d, consecutive failures: %"+fail_digits+"d, %"+fail_digits+"d.\r") % \
                             (len(unique_seeds), fails, recycle_fails))

        lightcone1.rerandomize_light_cone_solution(new_seed,
                                                   recycle=False)
        if new_recycle_seed is not None:
            lightcone1.rerandomize_light_cone_solution(new_recycle_seed,
                                                       recycle=True)

        # Compare with all other seeds.
        test_pass = True
        for unique_seed in unique_seeds:
            lightcone2.rerandomize_light_cone_solution(unique_seed['master'],
                                                       recycle=False)
            if unique_seed['recycle'] is not None:
                lightcone2.rerandomize_light_cone_solution(unique_seed['recycle'],
                                                           recycle=True)

            common = _compare_solutions(lightcone1.light_cone_solution,
                                        lightcone2.light_cone_solution)

            if (common > max_overlap):
                test_pass = False
                break
            else:
                max_common = max(max_common, common)

        if test_pass:
            unique_seeds.append({'master':new_seed,
                                 'recycle':new_recycle_seed})
            fails = 0
            recycle_fails = 0

        else:
            if recycle:
                recycle_fails += 1
            else:
                fails += 1

            if (recycle_fails >= failures):
                sys.stderr.write(("Unique solutions: %d, consecutive failures: %"+fail_digits+"d, %"+fail_digits+"d.\n") % \
                                 (len(unique_seeds), fails, recycle_fails))
                fails += 1
                mylog.info("Max recycled failures reached with master seed %d." % \
                           new_seed)
                master = None
            if (fails >= failures):
                sys.stderr.write(("Unique solutions: %d, consecutive failures: %"+fail_digits+"d, %"+fail_digits+"d.\n") % \
                                 (len(unique_seeds), fails, recycle_fails))
                mylog.error("Max consecutive failures reached.")
                break

    mylog.info("Created %d unique solutions." % len(unique_seeds))
    mylog.info("Maximum common volume is %.2e." % max_common)
    _write_seed_file(unique_seeds, filename)
    return unique_seeds

def _compare_solutions(solution1, solution2):
    "Calculate common volume between two light cone solutions."

    if (len(solution1) != len(solution2)):
        mylog.error("Cannot compare light cone solutions with unequal numbers of slices.")
        return -1

    my_volume = 0.0
    total_volume = 0.0

    # Check that solution volumes are the same.
    if((solution1[0]['box_depth_fraction'] *
        solution1[0]['box_width_fraction']**2) !=
       (solution2[0]['box_depth_fraction'] *
        solution2[0]['box_width_fraction']**2)):
        mylog.error("Light cone solutions do not have equal volumes, will use the smaller one.")

    for q in range(len(solution1)):
        cube1 = np.zeros(shape=(len(solution1[q]['projection_center']), 2))
        volume1 = 1.0
        for w in range(len(cube1)):
            if (w == solution1[q]['projection_axis']):
                width = solution1[q]['box_depth_fraction']
            else:
                width = solution1[q]['box_width_fraction']
            volume1 *= width
            cube1[w] = [solution1[q]['projection_center'][w] - 0.5 * width,
                        solution1[q]['projection_center'][w] + 0.5 * width]

        cube2 = np.zeros(shape=(len(solution2[q]['projection_center']), 2))
        volume2 = 1.0
        for w in range(len(cube2)):
            if (w == solution2[q]['projection_axis']):
                width = solution2[q]['box_depth_fraction']
            else:
                width = solution2[q]['box_width_fraction']
            volume2 *= width
            cube2[w] = [solution2[q]['projection_center'][w] - 0.5 * width,
                        solution2[q]['projection_center'][w] + 0.5 * width]

        total_volume += min(volume1, volume2)
        my_volume += common_volume(cube1, cube2,
                                   periodic=np.array([[0, 1],
                                                      [0, 1],
                                                      [0, 1]]))

    return (my_volume / total_volume)

def _read_seed_file(filename):
    "Read list of random seeds from a file."

    mylog.info("Reading random seed list from %s." % filename)

    seed_list = []

    lines = file(filename)
    for line in lines:
        if line[0] != '#':
            line = line.strip()
            onLine = line.split(',')
            if (len(onLine) == 1):
                seed_list.append({'master':onLine[0], 'recycle':None})
            else:
                seed_list.append({'master':onLine[0], 'recycle':onLine[1]})

    return seed_list

@parallel_root_only
def _write_seed_file(seed_list, filename):
    "Write list of random seeds to a file."

    mylog.info("Writing random seed list to %s." % filename)

    f = open(filename, 'w')
    for seed in seed_list:
        if seed['recycle'] is None:
            f.write("%s\n" % seed['master'])
        else:
            f.write("%s,%s\n" % (seed['master'], seed['recycle']))
    f.close()
