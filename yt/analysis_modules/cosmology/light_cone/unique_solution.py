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
import numpy as na
import random as rand
import sys

from yt.utilities.logger import ytLogger as mylog
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_root_only

from .light_cone import LightCone
from .common_n_volume import common_volume

def project_unique_light_cones(lightcone, seed_file, field, **kwargs):
    "Make light cone projections using list of random seeds in a file."

    seedList = _read_seed_file(seed_file)

    prefix = lightcone.output_prefix
    lightcone.calculate_light_cone_solution(seed=0)
    lastSeed = None

    for seed in seedList:
        if (seed['master'] != lastSeed):
            lightcone.rerandomize_light_cone_solution(seed['master'], recycle=False)
            lastSeed = seed['master']
        if (seed['recycle'] is not None):
            lightcone.rerandomize_light_cone_solution(seed['recycle'], recycle=True)

        lightcone.output_prefix = "%s_%s_%s" % (prefix, seed['master'], seed['recycle'])
        lightcone.project_light_cone(field, **kwargs)

def find_unique_solutions(lightcone1, solutions=100, seed=None, max_overlap=0.25, failures=10, 
                          recycle=True, filename='unique.dat'):
    "Find a set of random seeds that will give light cones will minimal volume overlap."

    lightcone2 = copy.deepcopy(lightcone1)
    lightcone1.calculate_light_cone_solution(seed=0)
    lightcone2.calculate_light_cone_solution(seed=0)

    uniqueSeeds = []
    if recycle:
        master = None
    newRecycleSeed = None
    fails = 0
    recycleFails = 0

    maxCommon = 0.0

    # Need to continuall save and reset the state of the random number generator
    # since it is being reset by the light cone generator.
    if seed is None:
        state = None
    else:
        rand.seed(seed)
        state = rand.getstate()

    failDigits = str(int(na.log10(failures))+1)

    while (len(uniqueSeeds) < solutions):
        # Create new random seed.
        if (recycle and master is not None):
            newSeed = master
            if state is not None: rand.setstate(state)
            newRecycleSeed = rand.randint(1, 1e9)
            state = rand.getstate()
        else:
            if state is not None: rand.setstate(state)
            newSeed = rand.randint(1, 1e9)
            state = rand.getstate()
            if recycle:
                master = newSeed
                recycleFails = 0
            newRecycleSeed = None

        sys.stderr.write(("Unique solutions: %d, consecutive failures: %"+failDigits+"d, %"+failDigits+"d.\r") % \
                             (len(uniqueSeeds), fails, recycleFails))

        lightcone1.rerandomize_light_cone_solution(newSeed, recycle=False)
        if newRecycleSeed is not None:
            lightcone1.rerandomize_light_cone_solution(newRecycleSeed, recycle=True)

        # Compare with all other seeds.
        testPass = True
        for uniqueSeed in uniqueSeeds:
            lightcone2.rerandomize_light_cone_solution(uniqueSeed['master'], recycle=False)
            if uniqueSeed['recycle'] is not None:
                lightcone2.rerandomize_light_cone_solution(uniqueSeed['recycle'], recycle=True)

            common = _compare_solutions(lightcone1.light_cone_solution, lightcone2.light_cone_solution)

            if (common > max_overlap):
                testPass = False
                break
            else:
                maxCommon = max(maxCommon, common)

        if testPass:
            uniqueSeeds.append({'master':newSeed, 'recycle':newRecycleSeed})
            fails = 0
            recycleFails = 0

        else:
            if recycle:
                recycleFails += 1
            else:
                fails += 1

            if (recycleFails >= failures):
                sys.stderr.write(("Unique solutions: %d, consecutive failures: %"+failDigits+"d, %"+failDigits+"d.\n") % \
                                     (len(uniqueSeeds), fails, recycleFails))
                fails += 1
                mylog.info("Max recycled failures reached with master seed %d." % newSeed)
                master = None
            if (fails >= failures):
                sys.stderr.write(("Unique solutions: %d, consecutive failures: %"+failDigits+"d, %"+failDigits+"d.\n") % \
                                     (len(uniqueSeeds), fails, recycleFails))
                mylog.error("Max consecutive failures reached.")
                break

    mylog.info("Created %d unique solutions." % len(uniqueSeeds))
    mylog.info("Maximum common volume is %.2e." % maxCommon)
    _write_seed_file(uniqueSeeds, filename)
    return uniqueSeeds

def _compare_solutions(solution1, solution2):
    "Calculate common volume between two light cone solutions."

    if (len(solution1) != len(solution2)):
        mylog.error("Cannot compare light cone solutions with unequal numbers of slices.")
        return -1

    commonVolume = 0.0
    totalVolume = 0.0

    # Check that solution volumes are the same.
    if((solution1[0]['box_depth_fraction'] * solution1[0]['box_width_fraction']**2) !=
       (solution2[0]['box_depth_fraction'] * solution2[0]['box_width_fraction']**2)):
        mylog.error("Light cone solutions do not have equal volumes, will use the smaller one.")

    for q in range(len(solution1)):
        cube1 = na.zeros(shape=(len(solution1[q]['projection_center']), 2))
        volume1 = 1.0
        for w in range(len(cube1)):
            if (w == solution1[q]['projection_axis']):
                width = solution1[q]['box_depth_fraction']
            else:
                width = solution1[q]['box_width_fraction']
            volume1 *= width
            cube1[w] = [solution1[q]['projection_center'][w] - 0.5 * width,
                        solution1[q]['projection_center'][w] + 0.5 * width]

        cube2 = na.zeros(shape=(len(solution2[q]['projection_center']), 2))
        volume2 = 1.0
        for w in range(len(cube2)):
            if (w == solution2[q]['projection_axis']):
                width = solution2[q]['box_depth_fraction']
            else:
                width = solution2[q]['box_width_fraction']
            volume2 *= width
            cube2[w] = [solution2[q]['projection_center'][w] - 0.5 * width,
                        solution2[q]['projection_center'][w] + 0.5 * width]

        totalVolume += min(volume1, volume2)
        commonVolume += common_volume(cube1, cube2, periodic=na.array([[0, 1], [0, 1], [0, 1]]))

    return (commonVolume/totalVolume)

def _read_seed_file(filename):
    "Read list of random seeds from a file."

    mylog.info("Reading random seed list from %s." % filename)

    seedList = []

    lines = file(filename)
    for line in lines:
        if line[0] != '#':
            line = line.strip()
            onLine = line.split(',')
            if (len(onLine) == 1):
                seedList.append({'master':onLine[0], 'recycle':None})
            else:
                seedList.append({'master':onLine[0], 'recycle':onLine[1]})

    return seedList

@parallel_root_only
def _write_seed_file(seedList, filename):
    "Write list of random seeds to a file."

    mylog.info("Writing random seed list to %s." % filename)

    f = open(filename, 'w')
    for seed in seedList:
        if seed['recycle'] is None:
            f.write("%s\n" % seed['master'])
        else:
            f.write("%s,%s\n" % (seed['master'], seed['recycle']))
    f.close()
