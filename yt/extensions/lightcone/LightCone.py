"""
LightCone class and member functions.

Author: Britton Smith <brittons@origins.colorado.edu>
Affiliation: CASA/University of CO, Boulder
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 Britton Smith.  All Rights Reserved.

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

from yt.extensions.lightcone import *
from yt.logger import lagosLogger as mylog
from yt.config import ytcfg
from yt.funcs import *
from Common_nVolume import *
from HaloMask import *
import copy
import os
import numpy as na
import random as rand

class LightCone(object):
    def __init__(self,EnzoParameterFile,LightConeParameterFile,verbose=True):
        self.verbose = verbose
        self.EnzoParameterFile = EnzoParameterFile
        self.LightConeParameterFile = LightConeParameterFile
        self.enzoParameters = {}
        self.lightConeParameters = {}
        self.redshiftOutputs = []
        self.timeOutputs = []
        self.allOutputs = []
        self.lightConeSolution = []
        self.masterSolution = [] # kept to compare with recycled solutions
        self.projectionStack = []
        self.projectionWeightFieldStack = []
        self.haloMask = []

        # Parameters for recycling light cone solutions.
        self.recycleSolution = False
        self.recycleRandomSeed = 0

        # Set some parameter defaults.
        self._SetParameterDefaults()

        # Read parameters.
        self._ReadLightConeParameterFile()
        self._ReadEnzoParameterFile()

        # Calculate number of pixels.
        self.pixels = int(self.lightConeParameters['FieldOfViewInArcMinutes'] * 60.0 / \
                              self.lightConeParameters['ImageResolutionInArcSeconds'])

        if ytcfg.getint("yt","__parallel_rank") == 0:
            # Create output directory.
            if (os.path.exists(self.lightConeParameters['OutputDir'])):
                if not(os.path.isdir(self.lightConeParameters['OutputDir'])):
                    mylog.error("Output directory exists, but is not a directory: %s." % self.lightConeParameters['OutputDir'])
                    self.lightConeParameters['OutputDir'] = './'
            else:
                os.mkdir(self.lightConeParameters['OutputDir'])

        # Calculate redshifts for dt data dumps.
        if (self.enzoParameters.has_key('dtDataDump')):
            self._CalculateTimeDumpRedshifts()

        # Combine all data dumps.
        self._CombineDataOutputs()

        # Calculate maximum delta z for each data dump.
        self._CalculateDeltaZMax()

        if not self.lightConeParameters['UseMinimumNumberOfProjections']:
            # Calculate minimum delta z for each data dump.
            self._CalculateDeltaZMin()

    def CalculateLightConeSolution(self,seed=None):
        "Create list of projections to be added together to make the light cone."

        # Make sure recycling flag is off.
        self.recycleSolution = False

        # Get rid of old halo mask, if one was there.
        self.haloMask = []

        if seed is not None:
            self.lightConeParameters['RandomSeed'] = int(seed)

        # Make light cone using minimum number of projections.
        if (self.lightConeParameters['UseMinimumNumberOfProjections']):

            z_Tolerance = 1e-4
            z = self.lightConeParameters['InitialRedshift']

            # fill redshift space with projections
            while ((z > self.lightConeParameters['FinalRedshift']) and 
                   (na.fabs(z - self.lightConeParameters['FinalRedshift']) > z_Tolerance)):
                # Sort data outputs by proximity to current redsfhit.
                self.allOutputs.sort(key=lambda obj:na.fabs(z - obj['redshift']))
                # For first data dump, choose closest to desired redshift.
                if (len(self.lightConeSolution) == 0):
                    self.lightConeSolution.append(self.allOutputs[0])
                # Start with data dump closest to desired redshift and move backward 
                # until one is within max delta z of last output in solution list.
                else:
                    output = self.allOutputs[0]
                    while (z > output['redshift']):
                        output = output['previous']
                        if (output is None):
                            if self.verbose: mylog.error("CalculateLightConeSolution: search for data output went off the end the stack.")
                            if self.verbose: mylog.error("Could not calculate light cone solution.")
                            return
                        if (output['redshift'] == self.lightConeSolution[-1]['redshift']):
                            if self.verbose: mylog.error("CalculateLightConeSolution: No data dump between z = %f and %f." % \
                                ((self.lightConeSolution[-1]['redshift'] - self.lightConeSolution[-1]['deltazMax']),
                                 self.lightConeSolution[-1]['redshift']))
                            if self.verbose: mylog.error("Could not calculate light cone solution.")
                            return
                    self.lightConeSolution.append(output)
                z = self.lightConeSolution[-1]['redshift'] - self.lightConeSolution[-1]['deltazMax']

        # Make light cone using maximum number of projections (minimum spacing).
        else:
            # Sort data outputs by proximity to current redsfhit.
            self.allOutputs.sort(key=lambda obj:na.fabs(self.lightConeParameters['InitialRedshift'] - obj['redshift']))
            # For first data dump, choose closest to desired redshift.
            self.lightConeSolution.append(self.allOutputs[0])

            nextOutput = self.lightConeSolution[-1]['next']
            while (nextOutput is not None):
                if (nextOutput['redshift'] <= self.lightConeParameters['FinalRedshift']):
                    break
                if ((self.lightConeSolution[-1]['redshift'] - nextOutput['redshift']) > self.lightConeSolution[-1]['deltazMin']):
                    self.lightConeSolution.append(nextOutput)
                nextOutput = nextOutput['next']

        if self.verbose: mylog.info("CalculateLightConeSolution: Used %d data dumps to get from z = %f to %f." % (len(self.lightConeSolution),
                                                                                            self.lightConeParameters['InitialRedshift'],
                                                                                            self.lightConeParameters['FinalRedshift']))

        # Calculate projection sizes, and get random projection axes and centers.
        rand.seed(self.lightConeParameters['RandomSeed'])
        co = lagos.Cosmology(HubbleConstantNow = (100.0 * self.enzoParameters['CosmologyHubbleConstantNow']),
                       OmegaMatterNow = self.enzoParameters['CosmologyOmegaMatterNow'],
                       OmegaLambdaNow = self.enzoParameters['CosmologyOmegaLambdaNow'])

        for q in range(len(self.lightConeSolution)):
            del self.lightConeSolution[q]['previous']
            del self.lightConeSolution[q]['next']
            if (q == len(self.lightConeSolution) - 1):
                z_next = self.lightConeParameters['FinalRedshift']
            else:
                z_next = self.lightConeSolution[q+1]['redshift']
            # Calculate fraction of box required for a depth of delta z
            self.lightConeSolution[q]['DepthBoxFraction'] = co.ComovingRadialDistance(z_next,self.lightConeSolution[q]['redshift']) * \
                self.enzoParameters['CosmologyHubbleConstantNow'] / self.enzoParameters['CosmologyComovingBoxSize']
            # Simple error check to make sure more than 100% of box depth is never required.
            if (self.lightConeSolution[q]['DepthBoxFraction'] > 1.0):
                if self.verbose: mylog.error("Warning: box fraction required to go from z = %f to %f is %f" % (self.lightConeSolution[q]['redshift'],z_next,
                                                                                                               self.lightConeSolution[q]['DepthBoxFraction']))
                if self.verbose: mylog.error("Full box delta z is %f, but it is %f to the next data dump." % (self.lightConeSolution[q]['deltazMax'],
                                                                                                              self.lightConeSolution[q]['redshift']-z_next))

            # Calculate fraction of box required for width corresponding to requested image size.
            scale = co.AngularScale_1arcsec_kpc(self.lightConeParameters['FinalRedshift'],self.lightConeSolution[q]['redshift'])
            size = self.lightConeParameters['FieldOfViewInArcMinutes'] * 60.0 * scale / 1000.0
            boxSizeProper = self.enzoParameters['CosmologyComovingBoxSize'] / (self.enzoParameters['CosmologyHubbleConstantNow'] * 
                                                                               (1.0 + self.lightConeSolution[q]['redshift']))
            self.lightConeSolution[q]['WidthBoxFraction'] = size / boxSizeProper

            # Get random projection axis and center.
            self.lightConeSolution[q]['ProjectionAxis'] = rand.randint(0,2)
            self.lightConeSolution[q]['ProjectionCenter'] = [rand.random(),rand.random(),rand.random()]

        # Store this as the master solution.
        self.masterSolution = [copy.deepcopy(q) for q in self.lightConeSolution]

        # Clear out some stuff.
        del co

    def GetHaloMask(self,HaloMaskParameterFile,mask_file=None,**kwargs):
        "Gets a halo mask from a file or makes a new one."

        # Check if file already exists.
        if (mask_file is not None) and os.path.exists(mask_file):
            input = h5.openFile(mask_file,'r')
            self.haloMask = input.root.haloMask.read()
            input.close()

        # Otherwise, make a halo mask.
        else:
            haloMaskCube = MakeLightConeHaloMask(self,HaloMaskParameterFile,mask_file=mask_file,**kwargs)
            # Collapse cube into final mask.
            self.haloMask = na.ones(shape=(self.pixels,self.pixels),dtype=bool)
            for mask in haloMaskCube:
                self.haloMask *= mask
            del haloMaskCube

    def ProjectLightCone(self,field,weight_field=None,apply_halo_mask=False,node=None,save_stack=True,save_slice_images=False,flatten_stack=False,**kwargs):
        "Create projections for light cone, then add them together."

        # Clear projection stack.
        self.projectionStack = []
        self.projectionWeightFieldStack = []
        if (self.lightConeSolution[-1].has_key('object')):
            del self.lightConeSolution[-1]['object']

        if not(self.lightConeParameters['OutputDir'].endswith("/")):
                 self.lightConeParameters['OutputDir'] += "/"

        for q,output in enumerate(self.lightConeSolution):
            if node is None:
                name = "%s%s_%04d_%04d" % (self.lightConeParameters['OutputDir'],self.lightConeParameters['OutputPrefix'],
                                           q,len(self.lightConeSolution))
            else:
                name = "%s%s_%s_%04d_%04d" % (self.lightConeParameters['OutputDir'],self.lightConeParameters['OutputPrefix'],
                                              node,q,len(self.lightConeSolution))
            output['object'] = lagos.EnzoStaticOutput(output['filename'])
            frb = LightConeProjection(output,field,self.pixels,weight_field=weight_field,
                                      save_image=save_slice_images,
                                      name=name,node=node,**kwargs)

            if ytcfg.getint("yt","__parallel_rank") == 0:
                if (weight_field is not None):
                    # Data come back normalized by the weight field.
                    # Undo that so it can be added up for the light cone.
                    self.projectionStack.append(frb[field]*frb['weight_field'])
                    self.projectionWeightFieldStack.append(frb['weight_field'])
                else:
                    self.projectionStack.append(frb[field])

                # Delete the frb.  This saves a decent amount of ram.
                if (q < len(self.lightConeSolution) - 1):
                    del frb

                # Flatten stack to save memory.
                if flatten_stack and (len(self.projectionStack) > 1):
                    self.projectionStack = [sum(self.projectionStack)]
                    if weight_field is not None:
                        self.projectionWeightFieldStack = [sum(self.projectionWeightFieldStack)]

            # Delete the plot collection now that the frb is deleted.
            del output['pc']

            # Unless this is the last slice, delete the dataset object.
            # The last one will be saved to make the plot collection.
            if (q < len(self.lightConeSolution) - 1):
                del output['object']

        if ytcfg.getint("yt","__parallel_rank") == 0:
            # Add up slices to make light cone projection.
            if (weight_field is None):
                lightConeProjection = sum(self.projectionStack)
            else:
                lightConeProjection = sum(self.projectionStack) / sum(self.projectionWeightFieldStack)

            if node is None:
                filename = "%s%s" % (self.lightConeParameters['OutputDir'],self.lightConeParameters['OutputPrefix'])
            else:
                filename = "%s%s_%s" % (self.lightConeParameters['OutputDir'],self.lightConeParameters['OutputPrefix'],node)

            # Save the last fixed resolution buffer for the plot collection, 
            # but replace the data with the full light cone projection data.
            frb.data[field] = lightConeProjection

            # Write stack to hdf5 file.
            if save_stack:
                self._SaveLightConeStack(field=field,weight_field=weight_field,filename=filename)

            # Apply halo mask.
            if apply_halo_mask:
                if len(self.haloMask) > 0:
                    mylog.info("Applying halo mask.")
                    frb.data[field] *= self.haloMask
                else:
                    mylog.error("No halo mask loaded, call GetHaloMask.")

            # Make a plot collection for the light cone projection.
            center = [0.5 * (self.lightConeSolution[-1]['object'].parameters['DomainLeftEdge'][w] + 
                             self.lightConeSolution[-1]['object'].parameters['DomainRightEdge'][w])
                      for w in range(self.lightConeSolution[-1]['object'].parameters['TopGridRank'])]
            pc = raven.PlotCollection(self.lightConeSolution[-1]['object'],center=center)
            pc.add_fixed_resolution_plot(frb,field)
            pc.save(filename)

            # Return the plot collection so the user can remake the plot if they want.
            return pc

    def RerandomizeLightConeSolution(self,newSeed,recycle=True):
        """
        When making a projection for a light cone, only randomizations along the line of sight make any 
        given projection unique, since the lateral shifting and tiling is done after the projection is made.
        Therefore, multiple light cones can be made from a single set of projections by introducing different 
        lateral random shifts and keeping all the original shifts along the line of sight.
        This routine will take in a new random seed and rerandomize the parts of the light cone that do not contribute 
        to creating a unique projection object.  Additionally, this routine is built such that if the same random 
        seed is given for the rerandomizing, the solution will be identical to the original.

        This routine has now been updated to be a general solution rescrambler.  If the keyword recycle is set to 
        True, then it will recycle.  Otherwise, it will create a completely new solution.
        """

        # Get rid of old halo mask, if one was there.
        self.haloMask = []

        # Clean pf objects out of light cone solution.
        for slice in self.lightConeSolution:
            if slice.has_key('object'):
                del slice['object']

        if recycle:
            if self.verbose: mylog.info("Recycling solution made with %s with new seed %s." % (self.lightConeParameters['RandomSeed'],
                                                                              newSeed))
            self.recycleRandomSeed = int(newSeed)
        else:
            if self.verbose: mylog.info("Creating new solution with random seed %s." % newSeed)
            self.lightConeParameters['RandomSeed'] = int(newSeed)
            self.recycleRandomSeed = 0

        self.recycleSolution = recycle

        # Keep track of fraction of volume in common between the original and recycled solution.
        commonVolume = 0.0
        totalVolume = 0.0

        # Seed random number generator with new seed.
        rand.seed(int(newSeed))

        for q,output in enumerate(self.lightConeSolution):
            # It is necessary to make the same number of calls to the random number generator
            # so the original solution willbe produced if the same seed is given.

            # Get random projection axis and center.
            # If recycling, axis will get thrown away since it is used in creating a unique projection object.
            newAxis = rand.randint(0,2)
            if recycle:
                output['ProjectionAxis'] = self.masterSolution[q]['ProjectionAxis']
            else:
                output['ProjectionAxis'] = newAxis

            newCenter = [rand.random(),rand.random(),rand.random()]

            # Make list of rectangle corners to calculate common volume.
            newCube = na.zeros(shape=(len(newCenter),2))
            oldCube = na.zeros(shape=(len(newCenter),2))
            for w in range(len(newCenter)):
                if (w == self.masterSolution[q]['ProjectionAxis']):
                    oldCube[w] = [self.masterSolution[q]['ProjectionCenter'][w] - 0.5 * self.masterSolution[q]['DepthBoxFraction'],
                                  self.masterSolution[q]['ProjectionCenter'][w] + 0.5 * self.masterSolution[q]['DepthBoxFraction']]
                else:
                    oldCube[w] = [self.masterSolution[q]['ProjectionCenter'][w] - 0.5 * self.masterSolution[q]['WidthBoxFraction'],
                                  self.masterSolution[q]['ProjectionCenter'][w] + 0.5 * self.masterSolution[q]['WidthBoxFraction']]

                if (w == output['ProjectionAxis']):
                    if recycle:
                        newCube[w] = oldCube[w]
                    else:
                        newCube[w] = [newCenter[w] - 0.5 * self.masterSolution[q]['DepthBoxFraction'],
                                      newCenter[w] + 0.5 * self.masterSolution[q]['DepthBoxFraction']]
                else:
                    newCube[w] = [newCenter[w] - 0.5 * self.masterSolution[q]['WidthBoxFraction'],
                                  newCenter[w] + 0.5 * self.masterSolution[q]['WidthBoxFraction']]

            commonVolume += commonNVolume(oldCube,newCube,periodic=na.array([[0,1],[0,1],[0,1]]))
            totalVolume += output['DepthBoxFraction'] * output['WidthBoxFraction']**2

            # Replace centers for every axis except the line of sight axis.
            for w in range(len(newCenter)):
                if not(recycle and (w == self.lightConeSolution[q]['ProjectionAxis'])):
                    self.lightConeSolution[q]['ProjectionCenter'][w] = newCenter[w]

        if recycle:
            if self.verbose: mylog.info("Fractional common volume between master and recycled solution is %.2e" % (commonVolume/totalVolume))
        else:
            if self.verbose: mylog.info("Fraction of total volume in common with old solution is %.2e." % (commonVolume/totalVolume))
            self.masterSolution = [copy.deepcopy(q) for q in self.lightConeSolution]

    def RestoreMasterSolution(self):
        "Reset the active light cone solution to the master solution."
        self.lightConeSolution = [copy.deepcopy(q) for q in self.masterSolution]

    @rootonly
    def SaveLightConeSolution(self,file="light_cone.dat"):
        "Write out a text file with information on light cone solution."

        if self.verbose: mylog.info("Saving light cone solution to %s." % file)

        f = open(file,'w')
        if self.recycleSolution:
            f.write("Recycled Solution\n")
            f.write("RecycleRandomSeed = %s\n" % self.recycleRandomSeed)
        else:
            f.write("Original Solution\n")
        f.write("EnzoParameterFile = %s\n" % self.EnzoParameterFile)
        f.write("LightConeParameterFile = %s\n" % self.LightConeParameterFile)
        f.write("\n")
        for parameter in self.lightConeParameters:
            f.write("%s = %s\n" % (parameter,str(self.lightConeParameters[parameter])))
        f.write("\n")
        for q,output in enumerate(self.lightConeSolution):
            f.write("Proj %04d, %s, z = %f, depth/box = %f, width/box = %f, axis = %d, center = %f, %f, %f\n" %
                    (q,output['filename'],output['redshift'],output['DepthBoxFraction'],output['WidthBoxFraction'],
                    output['ProjectionAxis'],output['ProjectionCenter'][0],output['ProjectionCenter'][1],output['ProjectionCenter'][2]))
        f.close()

    def _CalculateDeltaZMax(self):
        "Calculate delta z that corresponds to full box length going from z to (z - delta z)."
        co = lagos.Cosmology(HubbleConstantNow = (100.0 * self.enzoParameters['CosmologyHubbleConstantNow']),
                       OmegaMatterNow = self.enzoParameters['CosmologyOmegaMatterNow'],
                       OmegaLambdaNow = self.enzoParameters['CosmologyOmegaLambdaNow'])

        d_Tolerance = 1e-4
        max_Iterations = 100

        targetDistance = self.enzoParameters['CosmologyComovingBoxSize']

        for output in self.allOutputs:
            z = output['redshift']

            # Calculate delta z that corresponds to the length of the box at a given redshift.
            # Use Newton's method to calculate solution.
            z1 = z
            z2 = z1 - 0.1 # just an initial guess
            distance1 = 0.0
            iteration = 1

            # Convert comoving radial distance into Mpc / h, since that's how box size is stored.
            distance2 = co.ComovingRadialDistance(z2,z) * self.enzoParameters['CosmologyHubbleConstantNow']

            while ((na.fabs(distance2-targetDistance)/distance2) > d_Tolerance):
                m = (distance2 - distance1) / (z2 - z1)
                z1 = z2
                distance1 = distance2
                z2 = ((targetDistance - distance2) / m) + z2
                distance2 = co.ComovingRadialDistance(z2,z) * self.enzoParameters['CosmologyHubbleConstantNow']
                iteration += 1
                if (iteration > max_Iterations):
                    if self.verbose: mylog.error("CalculateDeltaZMax: Warning - max iterations exceeded for z = %f (delta z = %f)." % (z,na.fabs(z2-z)))
                    break
            output['deltazMax'] = na.fabs(z2-z)

        del co

    def _CalculateDeltaZMin(self):
        "Calculate delta z that corresponds to a single top grid pixel going from z to (z - delta z)."
        co = lagos.Cosmology(HubbleConstantNow = (100.0 * self.enzoParameters['CosmologyHubbleConstantNow']),
                       OmegaMatterNow = self.enzoParameters['CosmologyOmegaMatterNow'],
                       OmegaLambdaNow = self.enzoParameters['CosmologyOmegaLambdaNow'])

        d_Tolerance = 1e-4
        max_Iterations = 100

        targetDistance = self.enzoParameters['CosmologyComovingBoxSize'] / self.enzoParameters['TopGridDimensions'][0]

        for output in self.allOutputs:
            z = output['redshift']

            # Calculate delta z that corresponds to the length of a top grid pixel at a given redshift.
            # Use Newton's method to calculate solution.
            z1 = z
            z2 = z1 - 0.01 # just an initial guess
            distance1 = 0.0
            iteration = 1

            # Convert comoving radial distance into Mpc / h, since that's how box size is stored.
            distance2 = co.ComovingRadialDistance(z2,z) * self.enzoParameters['CosmologyHubbleConstantNow']

            while ((na.fabs(distance2 - targetDistance) / distance2) > d_Tolerance):
                m = (distance2 - distance1) / (z2 - z1)
                z1 = z2
                distance1 = distance2
                z2 = ((targetDistance - distance2) / m) + z2
                distance2 = co.ComovingRadialDistance(z2,z) * self.enzoParameters['CosmologyHubbleConstantNow']
                iteration += 1
                if (iteration > max_Iterations):
                    if self.verbose: mylog.error("CalculateDeltaZMax: Warning - max iterations exceeded for z = %f (delta z = %f)." % (z,na.fabs(z2-z)))
                    break
            output['deltazMin'] = na.fabs(z2-z)

        del co

    def _CalculateTimeDumpRedshifts(self):
        "Calculates redshift values for the dt data dumps."
        ec = lagos.EnzoCosmology(HubbleConstantNow = \
                               (100.0 * self.enzoParameters['CosmologyHubbleConstantNow']),
                           OmegaMatterNow = self.enzoParameters['CosmologyOmegaMatterNow'],
                           OmegaLambdaNow = self.enzoParameters['CosmologyOmegaLambdaNow'],
                           InitialRedshift = self.enzoParameters['CosmologyInitialRedshift'])

        z_Tolerance = 1e-3
        z = self.enzoParameters['CosmologyInitialRedshift']
        index = 0
        while ((z > self.enzoParameters['CosmologyFinalRedshift']) and 
               (na.fabs(z-self.enzoParameters['CosmologyFinalRedshift']) > z_Tolerance)):
            t = ec.InitialTime + (self.enzoParameters['dtDataDump'] * ec.TimeUnits * index)
            z = ec.ComputeRedshiftFromTime(t)
            filename = "%s/%s%04d/%s%04d" % (self.enzoParameters['GlobalDir'],
                                             self.enzoParameters['DataDumpDir'],
                                             index,
                                             self.enzoParameters['DataDumpDir'],
                                             index)
            self.timeOutputs.append({'index':index,'redshift':z,'filename':filename})
            index += 1

    def _CombineDataOutputs(self):
        "Combines redshift and time data into one sorted list."
        self.allOutputs = self.redshiftOutputs + self.timeOutputs
        self.allOutputs.sort(reverse=True,key=lambda obj:obj['redshift'])
        for q in range(len(self.allOutputs)):
            del self.allOutputs[q]['index']
            if (q == 0):
                self.allOutputs[q]['previous'] = None
                self.allOutputs[q]['next'] = self.allOutputs[q+1]
            elif (q == len(self.allOutputs) - 1):
                self.allOutputs[q]['previous'] = self.allOutputs[q-1]                
                self.allOutputs[q]['next'] = None
            else:
                self.allOutputs[q]['previous'] = self.allOutputs[q-1]
                self.allOutputs[q]['next'] = self.allOutputs[q+1]
        del self.redshiftOutputs
        del self.timeOutputs

    def _ReadEnzoParameterFile(self):
        "Reads an Enzo parameter file looking for cosmology and output parameters."
        lines = open(self.EnzoParameterFile).readlines()
        for line in lines:
            if line.find("#") >= 0: # Keep the commented lines
                line=line[:line.find("#")]
            line=line.strip()
            if len(line) < 2:
                continue
            try:
                param, vals = map(str,line.split("="))
                param = param.strip()
                vals = vals.strip()
            except ValueError:
                if self.verbose: mylog.error("Skipping line: %s" % line)
                continue
            if EnzoParameterDict.has_key(param):
                t = map(EnzoParameterDict[param], vals.split())
                if len(t) == 1:
                    self.enzoParameters[param] = t[0]
                else:
                    self.enzoParameters[param] = t
            elif param.startswith("CosmologyOutputRedshift["):
                index = param[param.find("[")+1:param.find("]")]
                self.redshiftOutputs.append({'index':int(index),
                                             'redshift':float(vals)})

        # Add filenames to redshift outputs.
        for output in self.redshiftOutputs:
            output["filename"] = "%s/%s%04d/%s%04d" % (self.enzoParameters['GlobalDir'],
                                                       self.enzoParameters['RedshiftDumpDir'],
                                                       output['index'],
                                                       self.enzoParameters['RedshiftDumpDir'],
                                                       output['index'])

    def _ReadLightConeParameterFile(self):
        "Reads a light cone parameter file."
        lines = open(self.LightConeParameterFile).readlines()
        for line in lines:
            if line.find("#") >= 0: # Keep the commented lines
                line=line[:line.find("#")]
            line=line.strip()
            if len(line) < 2:
                continue
            try:
                param, vals = map(str,line.split("="))
                param = param.strip()
                vals = vals.strip()
            except ValueError:
                if self.verbose: mylog.error("Skipping line: %s" % line)
                continue
            if LightConeParameterDict.has_key(param):
                t = map(LightConeParameterDict[param], vals.split())
                if len(t) == 1:
                    self.lightConeParameters[param] = t[0]
                else:
                    self.lightConeParameters[param] = t

    def _SaveLightConeStack(self,field=None,weight_field=None,filename=None):
        "Save the light cone projection stack as a 3d array in and hdf5 file."

        field_node = "%s_%s" % (field,weight_field)
        weight_field_node = "weight_field_%s" % weight_field

        import tables
        if (filename is None):
            filename = "%s/%s_data" % (self.lightConeParameters['OutputDir'],self.lightConeParameters['OutputPrefix'])
        if not(filename.endswith('.h5')):
               filename += ".h5"

        if (len(self.projectionStack) == 0):
            if self.verbose: mylog.error("SaveLightConeStack: no projection data loaded.")
            return

        if self.verbose: mylog.info("Writing light cone data to %s." % filename)

        output = tables.openFile(filename, "a")

        try:
            node_exists = output.isVisibleNode("/%s" % field_node)
        except tables.exceptions.NoSuchNodeError:
            node_exists = False

        if node_exists:
            mylog.error("Dataset, %s, already exists in %s, not saving." % (field_node,filename))
        else:
            mylog.info("Saving %s to %s." % (field_node, filename))
            self.projectionStack = na.array(self.projectionStack)
            output.createArray("/",field_node,self.projectionStack)

        if (len(self.projectionWeightFieldStack) > 0):
            try:
                node_exists = output.isVisibleNode("/%s" % weight_field_node)
            except tables.exceptions.NoSuchNodeError:
                node_exists = False

            if node_exists:
                mylog.error("Dataset, %s, already exists in %s, not saving." % (weight_field_node,filename))
            else:
                mylog.info("Saving %s to %s." % (weight_field_node, filename))
                self.projectionWeightFieldStack = na.array(self.projectionWeightFieldStack)
                output.createArray("/",weight_field_node,self.projectionWeightFieldStack)

        output.close()

    def _SetParameterDefaults(self):
        "Set some default parameters to avoid problems if they are not in the parameter file."
        self.enzoParameters['GlobalDir'] = ""
        self.enzoParameters['RedshiftDumpName'] = "RD"
        self.enzoParameters['RedshiftDumpDir'] = "RD"
        self.enzoParameters['DataDumpName'] = "DD"
        self.enzoParameters['DataDumpDir'] = "DD"
        self.lightConeParameters['UseMinimumNumberOfProjections'] = 1
        self.lightConeParameters['OutputDir'] = "./"
        self.lightConeParameters['OutputPrefix'] = "LightCone"

EnzoParameterDict = {"CosmologyCurrentRedshift": float,
                     "CosmologyComovingBoxSize": float,
                     "CosmologyOmegaMatterNow": float,
                     "CosmologyOmegaLambdaNow": float,
                     "CosmologyHubbleConstantNow": float,
                     "CosmologyInitialRedshift": float,
                     "CosmologyFinalRedshift": float,
                     "TopGridDimensions": float,
                     "dtDataDump": float,
                     "RedshiftDumpName": str,
                     "RedshiftDumpDir":  str,
                     "DataDumpName": str,
                     "DataDumpDir": str,
                     "GlobalDir" : str}

LightConeParameterDict = {"InitialRedshift": float,
                          "FinalRedshift": float,
                          "FieldOfViewInArcMinutes": float,
                          "ImageResolutionInArcSeconds": float,
                          "RandomSeed": int,
                          "UseMinimumNumberOfProjections": int,
                          "OutputDir": str,
                          "OutputPrefix": str}
