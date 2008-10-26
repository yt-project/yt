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

from yt.lagos.lightcone import *
import numpy as na
import random as rand
import tables as h5

class LightCone(object):
    def __init__(self,EnzoParameterFile,LightConeParameterFile):
        self.EnzoParameterFile = EnzoParameterFile
        self.LightConeParameterFile = LightConeParameterFile
        self.enzoParameters = {}
        self.lightConeParameters = {}
        self.redshiftOutputs = []
        self.timeOutputs = []
        self.allOutputs = []
        self.lightConeSolution = []
        self.projectionStack = []
        self.projectionWeightFieldStack = []

        # Set some parameter defaults.
        self._SetParameterDefaults()

        # Read parameters.
        self._ReadLightConeParameterFile()
        self._ReadEnzoParameterFile()

        # Calculate redshifts for dt data dumps.
        if (self.enzoParameters.has_key('dtDataDump')):
            self._CalculateTimeDumpRedshifts()

        # Combine all data dumps.
        self._CombineDataOutputs()

        # Calculate maximum delta z for each data dump.
        self._CalculateDeltaZMax()

    def CalculateLightConeSolution(self):
        "Create list of projections to be added together to make the light cone."

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
                # until one is within delta z of last output in solution list.
                else:
                    output = self.allOutputs[0]
                    while (z > output['redshift']):
                        output = output['previous']
                        if (output is None):
                            print "CalculateLightConeSolution: search for data output went off the end the stack."
                            print "Could not calculate light cone solution."
                            return
                        if (output['redshift'] == self.lightConeSolution[-1]['redshift']):
                            print "CalculateLightConeSolution: No data dump between z = %f and %f." % \
                                ((self.lightConeSolution[-1]['redshift'] - self.lightConeSolution[-1]['deltaz']),
                                 self.lightConeSolution[-1]['redshift'])
                            print "Could not calculate light cone solution."
                            return
                    self.lightConeSolution.append(output)
                z = self.lightConeSolution[-1]['redshift'] - self.lightConeSolution[-1]['deltaz']

        # Make light cone using maximum number of projections (minimum spacing).
        else:
            deltazMin = 0.01

            # Sort data outputs by proximity to current redsfhit.
            self.allOutputs.sort(key=lambda obj:na.fabs(self.lightConeParameters['InitialRedshift'] - obj['redshift']))
            # For first data dump, choose closest to desired redshift.
            self.lightConeSolution.append(self.allOutputs[0])

            nextOutput = self.lightConeSolution[-1]['next']
            while (nextOutput is not None):
                if (nextOutput['redshift'] <= self.lightConeParameters['FinalRedshift']):
                    break
                if ((self.lightConeSolution[-1]['redshift'] - nextOutput['redshift']) > deltazMin):
                    self.lightConeSolution.append(nextOutput)
                nextOutput = nextOutput['next']

        print "CalculateLightConeSolution: Used %d data dumps to get from z = %f to %f." % (len(self.lightConeSolution),
                                                                                            self.lightConeParameters['InitialRedshift'],
                                                                                            self.lightConeParameters['FinalRedshift'])

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
                print "Warning: box fraction required to go from z = %f to %f is %f" % (self.lightConeSolution[q]['redshift'],z_next,
                                                                self.lightConeSolution[q]['DepthBoxFraction'])
                print "Full box delta z is %f, but it is %f to the next data dump." % (self.lightConeSolution[q]['deltaz'],
                                                                                       self.lightConeSolution[q]['redshift']-z_next)

            # Calculate fraction of box required for width corresponding to requested image size.
            scale = co.AngularScale_1arcsec_kpc(self.lightConeParameters['FinalRedshift'],self.lightConeSolution[q]['redshift'])
            size = self.lightConeParameters['FieldOfViewInArcMinutes'] * 60.0 * scale / 1000.0
            boxSizeProper = self.enzoParameters['CosmologyComovingBoxSize'] / (self.enzoParameters['CosmologyHubbleConstantNow'] * 
                                                                               (1.0 + self.lightConeSolution[q]['redshift']))
            self.lightConeSolution[q]['WidthBoxFraction'] = size / boxSizeProper

            # Get random projection axis and center.
            self.lightConeSolution[q]['ProjectionAxis'] = rand.randint(0,2)
            self.lightConeSolution[q]['ProjectionCenter'] = [rand.random(),rand.random(),rand.random()]

        # Clear out some stuff.
        del self.allOutputs
        del self.redshiftOutputs
        del self.timeOutputs
        del co

    def ProjectLightCone(self,field,weight_field=None):
        "Create projections for light cone, then add them together."

        if not(self.lightConeParameters['OutputDir'].endswith("/")):
                 self.lightConeParameters['OutputDir'] += "/"

        pixels = int(self.lightConeParameters['FieldOfViewInArcMinutes'] * 60.0 / \
            self.lightConeParameters['ImageResolutionInArcSeconds'])

        for q,output in enumerate(self.lightConeSolution):
            name = "%s%s_%04d_%04d" % (self.lightConeParameters['OutputDir'],self.lightConeParameters['OutputPrefix'],
                                       q,len(self.lightConeSolution))
            output['object'] = EnzoStaticOutput(output['filename'])
            frb = LightConeProjection(output,field,pixels,weight_field=weight_field,
                                      save_image=self.lightConeParameters['SaveLightConeSlices'],name=name)
            if (weight_field is not None):
                # Data come back normalized by the weight field.
                # Undo that so it can be added up for the light cone.
                self.projectionStack.append(frb[field]*frb['weight_field'])
                self.projectionWeightFieldStack.append(frb['weight_field'])
            else:
                self.projectionStack.append(frb[field])

            # Unless this is the last slice, delete the dataset object.
            # The last one will be saved to make the plot collection.
            if (q < len(self.lightConeSolution) - 1):
                del output['object']

        # Add up slices to make light cone projection.
        if (weight_field is None):
            lightConeProjection = sum(self.projectionStack)
        else:
            lightConeProjection = sum(self.projectionStack) / sum(self.projectionWeightFieldStack)

        filename = "%s%s" % (self.lightConeParameters['OutputDir'],self.lightConeParameters['OutputPrefix'])

        # Save the last fixed resolution buffer for the plot collection, 
        # but replace the data with the full light cone projection data.
        frb.data[field] = lightConeProjection

        # Make a plot collection for the light cone projection.
        pc = raven.PlotCollection(self.lightConeSolution[-1]['object'],center=[0.5,0.5,0.5])
        pc.add_fixed_resolution_plot(frb,field)
        pc.save(filename)

        # Return the plot collection so the user can remake the plot if they want.
        return pc

    def SaveLightConeSolution(self,file="light_cone.dat"):
        "Write out a text file with information on light cone solution."

        print "Saving light cone solution to %s." % file

        f = open(file,'w')
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

    def SaveLightConeStack(self,filename=None):
        "Save the light cone projection stack as a 3d array in and hdf5 file."
        if (filename is None):
            filename = "%s%s_data" % (self.lightConeParameters['OutputDir'],self.lightConeParameters['OutputPrefix'])
        filename += ".h5"

        if (len(self.projectionStack) == 0):
            print "SaveLightConeStack: no projection data loaded."
            return

        print "Writing light cone data to %s." % filename

        output = h5.openFile(filename, "a")

        self.projectionStack = na.array(self.projectionStack)
        output.createArray("/", "data",self.projectionStack)

        if (len(self.projectionWeightFieldStack) > 0):
            self.projectionWeightFieldStack = na.array(self.projectionWeightFieldStack)
            output.createArray("/","weight",self.projectionWeightFieldStack)

        output.close()

    def _CalculateDeltaZMax(self):
        "Calculate delta z that corresponds to full box length going from z to (z - delta z)."
        co = lagos.Cosmology(HubbleConstantNow = (100.0 * self.enzoParameters['CosmologyHubbleConstantNow']),
                       OmegaMatterNow = self.enzoParameters['CosmologyOmegaMatterNow'],
                       OmegaLambdaNow = self.enzoParameters['CosmologyOmegaLambdaNow'])

        d_Tolerance = 1e-4
        max_Iterations = 100

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

            while ((na.fabs(distance2-self.enzoParameters['CosmologyComovingBoxSize'])/distance2) > d_Tolerance):
                m = (distance2 - distance1) / (z2 - z1)
                z1 = z2
                distance1 = distance2
                z2 = ((self.enzoParameters['CosmologyComovingBoxSize'] - distance2) / m) + z2
                distance2 = co.ComovingRadialDistance(z2,z) * self.enzoParameters['CosmologyHubbleConstantNow']
                iteration += 1
                if (iteration > max_Iterations):
                    print "CalculateDeltaZMax: Warning - max iterations exceeded for z = %f (delta z = %f)." % (z,na.fabs(z2-z))
                    break
            output['deltaz'] = na.fabs(z2-z)

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
                print "Skipping line: %s" % line
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
                print "Skipping line: %s" % line
                continue
            if LightConeParameterDict.has_key(param):
                t = map(LightConeParameterDict[param], vals.split())
                if len(t) == 1:
                    self.lightConeParameters[param] = t[0]
                else:
                    self.lightConeParameters[param] = t

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
                          "SaveLightConeSlices": int,
                          "OutputDir": str,
                          "OutputPrefix": str}
