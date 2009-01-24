"""
Functions to generate unique light cone solutions.

Author: Britton Smith <brittons@origins.colorado.edu>
Affiliation: CASA/University of CO, Boulder
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2009 Britton Smith.  All Rights Reserved.

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

from LightCone import *
from Common_nVolume import *
from yt.logger import lagosLogger as mylog
import numpy as na
import random as rand
import sys

def ProjectUniqueLightCones(EnzoParameterFile,LightConeParameterFile,seedFile,field,**kwargs):
    "Make light cone projections using list of random seeds in a file."

    seedList = _ReadSeedFile(seedFile)

    lc = LightCone(EnzoParameterFile,LightConeParameterFile)
    prefix = lc.lightConeParameters['OutputPrefix']
    lc.CalculateLightConeSolution(seed=0)
    lastSeed = None

    for seed in seedList:
        if (seed['master'] != lastSeed):
            lc.RerandomizeLightConeSolution(seed['master'],recycle=False)
            lastSeed = seed['master']
        if (seed['recycle'] is not None):
            lc.RerandomizeLightConeSolution(seed['recycle'],recycle=True)

        lc.lightConeParameters['OutputPrefix'] = "%s_%s_%s" % (prefix,seed['master'],seed['recycle'])
        lc.ProjectLightCone(field,**kwargs)

def FindUniqueSolutions(EnzoParameterFile,LightConeParameterFile,solutions=100,seed=None,
                        max_overlap=0.25,failures=10,recycle=True,filename='unique.dat'):
    "Find a set of random seeds that will give light cones will minimal volume overlap."

    solution1 = LightCone(EnzoParameterFile,LightConeParameterFile,verbose=False)
    solution2 = LightCone(EnzoParameterFile,LightConeParameterFile,verbose=False)
    solution1.CalculateLightConeSolution(seed=0)
    solution2.CalculateLightConeSolution(seed=0)

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
            newRecycleSeed = rand.randint(1,1e9)
            state = rand.getstate()
        else:
            if state is not None: rand.setstate(state)
            newSeed = rand.randint(1,1e9)
            state = rand.getstate()
            if recycle:
                master = newSeed
                recycleFails = 0
            newRecycleSeed = None

        sys.stderr.write(("Unique solutions: %d, consecutive failures: %"+failDigits+"d, %"+failDigits+"d.\r") % (len(uniqueSeeds),fails,recycleFails))

        solution1.RerandomizeLightConeSolution(newSeed,recycle=False)
        if newRecycleSeed is not None:
            solution1.RerandomizeLightConeSolution(newRecycleSeed,recycle=True)

        # Compare with all other seeds.
        testPass = True
        for uniqueSeed in uniqueSeeds:
            solution2.RerandomizeLightConeSolution(uniqueSeed['master'],recycle=False)
            if uniqueSeed['recycle'] is not None:
                solution2.RerandomizeLightConeSolution(uniqueSeed['recycle'],recycle=True)

            common = CompareSolutions(solution1.lightConeSolution,solution2.lightConeSolution)

            if (common > max_overlap):
                testPass = False
                break
            else:
                maxCommon = max(maxCommon,common)

        if testPass:
            uniqueSeeds.append({'master':newSeed,'recycle':newRecycleSeed})
            fails = 0
            recycleFails = 0

        else:
            if recycle:
                recycleFails += 1
            else:
                fails += 1

            if (recycleFails >= failures):
                sys.stderr.write(("Unique solutions: %d, consecutive failures: %"+failDigits+"d, %"+failDigits+"d.\n") % (len(uniqueSeeds),fails,recycleFails))
                fails += 1
                mylog.info("Max recycled failures reached with master seed %d." % newSeed)
                master = None
            if (fails >= failures):
                sys.stderr.write(("Unique solutions: %d, consecutive failures: %"+failDigits+"d, %"+failDigits+"d.\n") % (len(uniqueSeeds),fails,recycleFails))
                mylog.error("Max consecutive failures reached.")
                break

    mylog.info("Created %d unique solutions." % len(uniqueSeeds))
    mylog.info("Maximum common volume is %.2e." % maxCommon)
    _WriteSeedFile(uniqueSeeds,filename)
    return uniqueSeeds

def CompareSolutions(solution1,solution2):
    "Calculate common volume between two light cone solutions."

    if (len(solution1) != len(solution2)):
        mylog.error("Cannot compare light cone solutions with unequal numbers of slices.")
        return -1

    commonVolume = 0.0
    totalVolume = 0.0

    # Check that solution volumes are the same.
    if((solution1[0]['DepthBoxFraction'] * solution1[0]['WidthBoxFraction']**2) !=
       (solution2[0]['DepthBoxFraction'] * solution2[0]['WidthBoxFraction']**2)):
        mylog.error("Light cone solutions do not have equal volumes, will use the smaller one.")

    for q in range(len(solution1)):
        cube1 = na.zeros(shape=(len(solution1[q]['ProjectionCenter']),2))
        volume1 = 1.0
        for w in range(len(cube1)):
            if (w == solution1[q]['ProjectionAxis']):
                width = solution1[q]['DepthBoxFraction']
            else:
                width = solution1[q]['WidthBoxFraction']
            volume1 *= width
            cube1[w] = [solution1[q]['ProjectionCenter'][w] - 0.5 * width,
                        solution1[q]['ProjectionCenter'][w] + 0.5 * width]

        cube2 = na.zeros(shape=(len(solution2[q]['ProjectionCenter']),2))
        volume2 = 1.0
        for w in range(len(cube2)):
            if (w == solution2[q]['ProjectionAxis']):
                width = solution2[q]['DepthBoxFraction']
            else:
                width = solution2[q]['WidthBoxFraction']
            volume2 *= width
            cube2[w] = [solution2[q]['ProjectionCenter'][w] - 0.5 * width,
                        solution2[q]['ProjectionCenter'][w] + 0.5 * width]

        totalVolume += min(volume1,volume2)
        commonVolume += commonNVolume(cube1,cube2,periodic=na.array([[0,1],[0,1],[0,1]]))

    return (commonVolume/totalVolume)

def _ReadSeedFile(filename):
    "Read list of random seeds from a file."

    mylog.info("Reading random seed list from %s." % filename)

    seedList = []

    lines = file(filename)
    for line in lines:
        line = line.strip()
        onLine = line.split(',')
        if (len(onLine) == 1):
            seedList.append({'master':onLine[0], 'recycle':None})
        else:
            seedList.append({'master':onLine[0], 'recycle':onLine[1]})

    return seedList

def _WriteSeedFile(seedList,filename):
    "Write list of random seeds to a file."

    mylog.info("Writing random seed list to %s." % filename)

    f = open(filename,'w')
    for seed in seedList:
        if seed['recycle'] is None:
            f.write("%s\n" % seed['master'])
        else:
            f.write("%s,%s\n" % (seed['master'],seed['recycle']))
    f.close()
