"""
HaloProfiler class and member functions.

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

import yt.lagos as lagos
import yt.lagos.hop as hop
from yt.logger import lagosLogger as mylog
import yt.raven as raven
import numpy as na
import os
import tables as h5

PROFILE_RADIUS_THRESHOLD = 2

class HaloProfiler(object):
    def __init__(self,dataset,HaloProfilerParameterFile,halos='multiple',radius=0.1,hop_style='new'):
        self.dataset = dataset
        self.HaloProfilerParameterFile = HaloProfilerParameterFile
        self.haloProfilerParameters = {}
        self.profileFields = {}
        self.projectionFields = {}
        self.hopHalos = []
        self.virialQuantities = []

        # Set option to get halos from hop or single halo at density maximum.
        # multiple: get halos from hop
        # single: get single halo from density maximum
        self.halos = halos
        if not(self.halos is 'multiple' or self.halos is 'single'):
            mylog.error("Keyword, halos, must be either 'single' or 'multiple'.")
            return None

        # Set hop file style.
        # old: enzo_hop output.
        # new: yt hop output.
        self.hop_style = hop_style
        if not(self.hop_style is 'old' or self.hop_style is 'new'):
            mylog.error("Keyword, hop_style, must be either 'old' or 'new'.")
            return None

        if self.halos is 'single' or hop_style is 'old':
            self.haloRadius = radius

        # Set some parameter defaults.
        self._SetParameterDefaults()

        # Read parameter file.
        self._ReadHaloProfilerParameterFile()

        # Check validity for VelocityCenter parameter which toggles how the 
        # velocity is zeroed out for radial velocity profiles.
        if self.haloProfilerParameters['VelocityCenter'][0] == 'bulk':
            if self.haloProfilerParameters['VelocityCenter'][1] == 'halo' and \
                    self.halos is 'single':
                mylog.error("Parameter, VelocityCenter, must be set to 'bulk sphere' or 'max <field>' with halos flag set to 'single'.")
                return None
            if self.haloProfilerParameters['VelocityCenter'][1] == 'halo' and \
                    self.hop_style is 'old':
                mylog.error("Parameter, VelocityCenter, must be 'bulk sphere' for old style hop output files.")
                return None
            if not(self.haloProfilerParameters['VelocityCenter'][1] == 'halo' or 
                   self.haloProfilerParameters['VelocityCenter'][1] == 'sphere'):
                mylog.error("Second value of VelocityCenter must be either 'halo' or 'sphere' if first value is 'bulk'.")
                return None
        elif self.haloProfilerParameters['VelocityCenter'][0] == 'max':
            if self.halos is 'multiple':
                mylog.error("Getting velocity center from a max field value only works with halos='single'.")
                return None
        else:
            mylog.error("First value of parameter, VelocityCenter, must be either 'bulk' or 'max'.")
            return None

        # Create dataset object.
        self.pf = lagos.EnzoStaticOutput(self.dataset)

    def makeProfiles(self):
        "Make radial profiles for all halos on the list."

        # Get halo(s).
        if self.halos is 'single':
            v, center = self.pf.h.find_max('Density')
            singleHalo = {}
            singleHalo['center'] = center
            singleHalo['r_max'] = self.haloRadius * self.pf.units['mpc']
            self.hopHalos.append(singleHalo)
        elif self.halos is 'multiple':
            # Get hop data.
            self._LoadHopData()
        else:
            mylog.error("I don't know whether to get halos from hop or from density maximum.  This should not have happened.")
            return

        # Add profile fields necessary for calculating virial quantities.
        self._CheckForNeededProfileFields()

        outputDir = "%s/%s" % (self.pf.fullpath,self.haloProfilerParameters['ProfileOutputDir'])

        if (os.path.exists(outputDir)):
            if not(os.path.isdir(outputDir)):
                mylog.error("Output directory exists, but is not a directory: %s." % outputDir)
                return
        else:
            os.mkdir(outputDir)

        pbar = lagos.get_pbar("Profiling halos ", len(self.hopHalos))
        for q,halo in enumerate(self.hopHalos):
            filename = "%s/Halo_%04d_profile.dat" % (outputDir,q)

            r_min = 2*self.pf.h.get_smallest_dx() * self.pf['mpc']
            if (halo['r_max'] / r_min < PROFILE_RADIUS_THRESHOLD):
                mylog.error("Skipping halo with r_max / r_min = %f." % (halo['r_max']/r_min))
                continue

            sphere = self.pf.h.sphere(halo['center'],halo['r_max']/self.pf.units['mpc'])

            # Set velocity to zero out radial velocity profiles.
            if self.haloProfilerParameters['VelocityCenter'][0] == 'bulk':
                if self.haloProfilerParameters['VelocityCenter'][1] == 'halo':
                    sphere.set_field_parameter('bulk_velocity',halo['velocity'])
                elif self.haloProfilerParameters['VelocityCenter'][1] == 'sphere':
                    sphere.set_field_parameter('bulk_velocity',sphere.quantities['BulkVelocity']())
                else:
                    mylog.error("Invalid parameter: VelocityCenter.")
            elif self.haloProfilerParameters['VelocityCenter'][0] == 'max':
                max_grid,max_cell,max_value,max_location = self.pf.h.find_max_cell_location(self.haloProfilerParameters['VelocityCenter'][1])
                sphere.set_field_parameter('bulk_velocity',[max_grid['x-velocity'][max_cell],
                                                            max_grid['y-velocity'][max_cell],
                                                            max_grid['z-velocity'][max_cell]])

            profile = lagos.BinnedProfile1D(sphere,self.haloProfilerParameters['n_bins'],"RadiusMpc",
                                            r_min,halo['r_max'],
                                            log_space=True, lazy_reader=True)
            for field in self.profileFields.keys():
                profile.add_fields(field,weight=self.profileFields[field][0],
                                   accumulation=self.profileFields[field][1])

            self._AddActualOverdensity(profile)

            virial = self._CalculateVirialQuantities(profile)
            virial['center'] = self.hopHalos[q]['center']

            if (virial['TotalMassMsun'] < self.haloProfilerParameters['VirialMassCutoff']):
                self.virialQuantities.append(None)
            else:
                self.virialQuantities.append(virial)
                profile.write_out(filename, format='%0.6e')
            del profile

            # Temporary solution to memory leak.
            for g in pf.h.grids:
                g.clear_data()
            sphere.clear_data()

            del sphere
            pbar.update(q)

        pbar.finish()
        self._WriteVirialQuantities()

    def makeProjections(self,save_images=True,save_cube=True,**kwargs):
        "Make projections of all halos using specified fields."

        # Get virial quantities.
        self._LoadVirialData()

        # Set resolution for fixed resolution output.
        if save_cube:
            if (str(self.haloProfilerParameters['ProjectAtLevel']).find('max') >= 0):
                proj_level = self.pf.h.maxLevel
            else:
                proj_level = int(self.haloProfilerParameters['ProjectAtLevel'])
            proj_dx = self.pf.units['mpc'] / self.pf.parameters['TopGridDimensions'][0] / \
                (self.pf.parameters['RefineBy']**proj_level)
            projectionResolution = int(self.haloProfilerParameters['ProjectionWidth'] / proj_dx)

        outputDir = "%s/%s" % (self.pf.fullpath,self.haloProfilerParameters['ProjectionOutputDir'])

        if (os.path.exists(outputDir)):
            if not(os.path.isdir(outputDir)):
                mylog.error("Output directory exists, but is not a directory: %s." % outputDir)
                return
        else:
            os.mkdir(outputDir)

        center = [0.5 * (self.pf.parameters['DomainLeftEdge'][w] + self.pf.parameters['DomainRightEdge'][w])
                  for w in range(self.pf.parameters['TopGridRank'])]

        # Create a plot collection.
        pc = raven.PlotCollection(self.pf,center=center)

        for q,halo in enumerate(self.virialQuantities):
            if halo is None:
                continue
            # Check if region will overlap domain edge.
            # Using non-periodic regions is faster than using periodic ones.
            leftEdge = [(halo['center'][w] - 0.5 * self.haloProfilerParameters['ProjectionWidth']/self.pf.units['mpc'])
                        for w in range(len(halo['center']))]
            rightEdge = [(halo['center'][w] + 0.5 * self.haloProfilerParameters['ProjectionWidth']/self.pf.units['mpc'])
                         for w in range(len(halo['center']))]

            mylog.info("Projecting halo %04d in region: [%f, %f, %f] to [%f, %f, %f]." %
                       (q,leftEdge[0],leftEdge[1],leftEdge[2],rightEdge[0],rightEdge[1],rightEdge[2]))

            need_per = False
            for w in range(len(halo['center'])):
                if ((leftEdge[w] < self.pf.parameters['DomainLeftEdge'][w]) or
                    (rightEdge[w] > self.pf.parameters['DomainRightEdge'][w])):
                    need_per = True
                    break

            if need_per:
                region = self.pf.h.periodic_region(halo['center'],leftEdge,rightEdge)
            else:
                region = self.pf.h.region(halo['center'],leftEdge,rightEdge)

            # Make projections.
            for w in range(self.pf.parameters['TopGridRank']):
                # YT projections do not follow the right-hand rule.
                coords = range(3)
                del coords[w]
                x_axis = coords[0]
                y_axis = coords[1]

                for field in self.projectionFields.keys():
                    pc.add_projection(field,w,weight_field=self.projectionFields[field],source=region,**kwargs)

                # Set x and y limits, shift image if it overlaps domain boundary.
                if need_per:
                    ShiftProjections(self.pf,pc,halo['center'],center,w)
                    # Projection has now been shifted to center of box.
                    proj_left = [center[x_axis]-0.5 * self.haloProfilerParameters['ProjectionWidth']/self.pf.units['mpc'],
                                 center[y_axis]-0.5 * self.haloProfilerParameters['ProjectionWidth']/self.pf.units['mpc']]
                    proj_right = [center[x_axis]+0.5 * self.haloProfilerParameters['ProjectionWidth']/self.pf.units['mpc'],
                                  center[y_axis]+0.5 * self.haloProfilerParameters['ProjectionWidth']/self.pf.units['mpc']]
                else:
                    proj_left = [leftEdge[x_axis],leftEdge[y_axis]]
                    proj_right = [rightEdge[x_axis],rightEdge[y_axis]]

                pc.set_xlim(proj_left[0],proj_right[0])
                pc.set_ylim(proj_left[1],proj_right[1])

                # Save projection data to hdf5 file.
                if save_cube:
                    axes = ['x','y','z']
                    dataFilename = "%s/Halo_%04d_%s_data.h5" % (outputDir,q,axes[w])
                    mylog.info("Saving projection data to %s." % dataFilename)

                    output = h5.openFile(dataFilename, "a")
                    # Create fixed resolution buffer for each projection and write them out.
                    for e,field in enumerate(self.projectionFields.keys()):
                        frb = raven.FixedResolutionBuffer(pc.plots[e].data,(proj_left[0],proj_right[0],proj_left[1],proj_right[1]),
                                                          (projectionResolution,projectionResolution),
                                                          antialias=False)
                        output.createArray("/",field,frb[field])
                    output.close()

                if save_images:
                    pc.save("%s/Halo_%04d" % (outputDir,q))

                pc.clear_plots()

            del region

        del pc

    def _WriteVirialQuantities(self):
        "Write out file with halo centers and virial masses and radii."
        filename = "%s/%s" % (self.pf.fullpath,self.haloProfilerParameters['VirialQuantitiesOutputFile'])
        mylog.info("Writing virial quantities to %s." % filename)
        file = open(filename,'w')
        file.write("#Index\tx\ty\tz\tMass [Msolar]\tRadius [Mpc]\n")
        for q in range(len(self.virialQuantities)):
            if (self.virialQuantities[q] is not None):
                file.write("%04d %.10f %.10f %.10f %.6e %.6e\n" % (q,self.hopHalos[q]['center'][0],self.hopHalos[q]['center'][1],
                                                                   self.hopHalos[q]['center'][2],
                                                                   self.virialQuantities[q]['TotalMassMsun'],
                                                                   self.virialQuantities[q]['RadiusMpc']))
        file.close()

    def _AddActualOverdensity(self,profile):
        "Calculate overdensity from TotalMassMsun and CellVolume fields."

        rho_crit_now = 1.8788e-29 # g cm^-3
        Msun2g = 1.989e33
        rho_crit = rho_crit_now * ((1 + self.pf['CosmologyCurrentRedshift'])**3.0)

        profile['ActualOverdensity'] = (Msun2g * profile['TotalMassMsun']) / \
            profile['CellVolume'] / rho_crit

    def _CalculateVirialQuantities(self,profile,overdensity_field='ActualOverdensity'):
        "Calculate virial radius and virial mass from radial profiles."

        fields = ['TotalMassMsun','RadiusMpc']
        overDensity = []
        temp_profile = {}
        for field in fields:
            temp_profile[field] = []

        for q in range(len(profile[overdensity_field])):
            good = True
            if (profile[overdensity_field][q] != profile[overdensity_field][q]):
                good = False
                continue
            for field in fields:
                if (profile[field][q] != profile[field][q]):
                    good = False
                    break
            if good:
                overDensity.append(profile[overdensity_field][q])
                for field in fields:
                    temp_profile[field].append(profile[field][q])

        virial = {}
        if (len(overDensity) < 2):
            mylog.error("Skipping halo with no valid points in profile.")
            for field in fields:
                virial[field] = 0.0
            return virial

        if (overDensity[1] <= self.haloProfilerParameters['VirialOverdensity']):
            index = 0
        elif (overDensity[-1] >= self.haloProfilerParameters['VirialOverdensity']):
            index = -2
        else:
            for q in (na.array(range(len(overDensity)-2))+2):
                if (overDensity[q] < self.haloProfilerParameters['VirialOverdensity']):
                    index = q - 1
                    break

        for field in fields:
            slope = (temp_profile[field][index+1] - temp_profile[field][index]) / \
                (overDensity[index+1] - overDensity[index])
            value = slope * (self.haloProfilerParameters['VirialOverdensity'] - overDensity[index]) + \
                temp_profile[field][index]
            virial[field] = value

        return virial

    def _CheckForNeededProfileFields(self):
        "Make sure CellVolume and TotalMass fields are added so virial quantities can be calculated."
        if not(self.profileFields.has_key('CellVolume')):
            mylog.info("Adding CellVolume field to so virial quantities can be calculated")
            self.profileFields['CellVolume'] = [None,True]
        if not(self.profileFields.has_key('TotalMassMsun')):
            mylog.info("Adding TotalMassMsun field to so virial quantities can be calculated")
            self.profileFields['TotalMassMsun'] = [None,True]

    def _ReadHopFile(self,hopFile):
        "Read hop file to get halo information."
        mylog.info("Reading halo information from %s." % hopFile)
        hopLines = file(hopFile)

        for line in hopLines:
            line = line.strip()
            if not(line.startswith('#')):
                onLine = line.split()
                mass = float(onLine[1])
                if (mass >= self.haloProfilerParameters['VirialMassCutoff']):
                    center = [float(onLine[7]),float(onLine[8]),float(onLine[9])]
                    velocity = [float(onLine[10]),float(onLine[11]),float(onLine[12])]
                    r_max = float(onLine[13]) * self.pf.units['mpc']
                    halo = {'center': center, 'r_max': r_max, 'velocity': velocity}
                    self.hopHalos.append(halo)

        mylog.info("Loaded %d halos with total dark matter mass af at least %e Msolar." % 
                   (len(self.hopHalos),self.haloProfilerParameters['VirialMassCutoff']))

    def _ReadOldHopFile(self,hopFile):
        "Read old style hop file made by enzo_hop."
        mylog.info("Reading halo information from old style hop file %s." % hopFile)
        hopLines = file(hopFile)

        for line in hopLines:
            line = line.strip()
            if not(line.startswith('#')):
                onLine = line.split()
                center = [float(onLine[4]),float(onLine[5]),float(onLine[6])]
                r_max = self.haloRadius * self.pf.units['mpc']
                halo = {'center': center, 'r_max': r_max}
                self.hopHalos.append(halo)

    def _RunHop(self,hopFile):
        "Run hop to get halos."
        full_box = self.pf.h.region(0.5 *(self.pf.parameters['DomainRightEdge'] -
                                          self.pf.parameters['DomainLeftEdge']),
                                    self.pf.parameters['DomainLeftEdge'],
                                    self.pf.parameters['DomainRightEdge'])

        hop_results = hop.HopList(full_box, 80.0)
        hop_results.write_out(hopFile)

        del full_box
        del hop_results

    def _LoadHopData(self):
        "Read hop output file or run hop if it doesn't exist."

        # Don't run if hop data already loaded.
        if self.hopHalos:
            return

        hopFile = "%s/%s" % (self.pf.fullpath,
                             self.haloProfilerParameters['HopOutputFile'])

        if not(os.path.exists(hopFile)):
            mylog.info("Hop file not found, running hop to get halos.")
            self._RunHop(hopFile)

        if self.hop_style is 'new':
            self._ReadHopFile(hopFile)
        else:
            self._ReadOldHopFile(hopFile)

    def _LoadVirialData(self):
        "Read virial quantities data or ask for profiles to be run if it doesn't exist."

        if self.virialQuantities:
            return

        virialFile = "%s/%s" % (self.pf.fullpath,
                                self.haloProfilerParameters['VirialQuantitiesOutputFile'])
        if not(os.path.exists(virialFile)):
            mylog.info("Virial quantities file not found.  Making profiles to calculate virial quantities.")
            self.makeProfiles()

        mylog.info("Reading virial quantities from %s." % virialFile)
        virialLines = file(virialFile)

        halos = 0
        for line in virialLines:
            line = line.strip()
            if not(line.startswith('#')):
                onLine = line.split()
                index = int(onLine[0])
                virial = {}
                virial['center'] =  [float(onLine[1]),float(onLine[2]),float(onLine[3])]
                virial['TotalMassMsun'] = float(onLine[4])
                virial['RadiusMpc'] = float(onLine[5])
                if (virial['TotalMassMsun'] >= self.haloProfilerParameters['VirialMassCutoff']):
                    for q in range(index - len(self.virialQuantities)):
                        self.virialQuantities.append(None)
                    self.virialQuantities.append(virial)
                    halos += 1

        mylog.info("Loaded virial quantities for %d halos." % halos)

    def _ReadHaloProfilerParameterFile(self):
        "Read a halo profiler parameter file."
        lines = open(self.HaloProfilerParameterFile).readlines()
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
                mylog.error("Skipping line: %s" % line)
                continue
            if haloProfilerParameterDict.has_key(param):
                t = map(haloProfilerParameterDict[param], vals.split())
                if len(t) == 1:
                    self.haloProfilerParameters[param] = t[0]
                else:
                    self.haloProfilerParameters[param] = t
            elif param.startswith("Profile["):
                field = param[param.find("[")+1:param.find("]")]
                field_vals = vals.split(',')
                for val in field_vals:
                    val.strip()
                if (field_vals[0].find('None') >= 0):
                    field_vals[0] = None
                if (field_vals[1].find('True') >= 0):
                    field_vals[1] = True
                else:
                    field_vals[1] = False
                self.profileFields[field] = field_vals
            elif param.startswith("Projection["):
                field = param[param.find("[")+1:param.find("]")]
                if (vals.find('None') >= 0):
                    vals = None
                self.projectionFields[field] = vals

    def _SetParameterDefaults(self):
        "Set some default parameters."
        self.haloProfilerParameters['VirialMassCutoff'] = 1e13 # M_solar
        self.haloProfilerParameters['VirialOverdensity'] = 200
        self.haloProfilerParameters['VelocityCenter'] = ['bulk','halo']
        self.haloProfilerParameters['n_bins'] = 50
        self.haloProfilerParameters['ProfileOutputDir'] = 'radial_profiles'
        self.haloProfilerParameters['ProjectionOutputDir'] = 'projections'
        self.haloProfilerParameters['HopOutputFile'] = "HopAnalysis.out"
        self.haloProfilerParameters['VirialQuantitiesOutputFile'] = "VirialQuantities.out"
        self.haloProfilerParameters['ProjectionWidth'] = 4.0 # Mpc
        self.haloProfilerParameters['ProjectAtLevel'] = 'max'

haloProfilerParameterDict = {"ProfileOutputDir": str,
                             "VirialQuantitiesOutputFile": str,
                             "HopOutputFile": str,
                             "VirialMassCutoff": float,
                             "VirialOverdensity": float,
                             "VelocityCenter": str,
                             "n_bins": int,
                             "ProjectionOutputDir": str,
                             "ProjectionWidth": float,
                             "ProjectAtLevel": str}

def ShiftProjections(pf,pc,oldCenter,newCenter,axis):
    """
    Shift projection data around.
    This is necessary when projecting a preiodic region.
    """
    offset = [newCenter[q]-oldCenter[q] for q in range(len(oldCenter))]
    width = [pf.parameters['DomainRightEdge'][q]-pf.parameters['DomainLeftEdge'][q] for q in range(len(oldCenter))]

    del offset[axis]
    del width[axis]

    for plot in pc.plots:
        # Get name of data field.
        other_fields = {'px':True,'py':True,'pdx':True,'pdy':True,'weight_field':True}
        for pfield in plot.data.data.keys():
            if not(other_fields.has_key(pfield)):
                field = pfield
                break

        # Shift x and y positions.
        plot['px'] += offset[0]
        plot['py'] += offset[1]

        # Wrap off-edge cells back around to other side (periodic boundary conditions).
        plot['px'][plot['px'] < 0] += width[0]
        plot['py'][plot['py'] < 0] += width[1]
        plot['px'][plot['px'] > width[0]] -= width[0]
        plot['py'][plot['py'] > width[1]] -= width[1]

        # After shifting, some cells have fractional coverage on both sides of the box.
        # Find those cells and make copies to be placed on the other side.

        # Cells hanging off the right edge.
        add_x_px = plot['px'][plot['px'] + 0.5 * plot['pdx'] > width[0]]
        add_x_px -= width[0]
        add_x_py = plot['py'][plot['px'] + 0.5 * plot['pdx'] > width[0]]
        add_x_pdx = plot['pdx'][plot['px'] + 0.5 * plot['pdx'] > width[0]]
        add_x_pdy = plot['pdy'][plot['px'] + 0.5 * plot['pdx'] > width[0]]
        add_x_field = plot[field][plot['px'] + 0.5 * plot['pdx'] > width[0]]
        add_x_weight_field = plot['weight_field'][plot['px'] + 0.5 * plot['pdx'] > width[0]]

        # Cells hanging off the left edge.
        add2_x_px = plot['px'][plot['px'] - 0.5 * plot['pdx'] < 0]
        add2_x_px += width[0]
        add2_x_py = plot['py'][plot['px'] - 0.5 * plot['pdx'] < 0]
        add2_x_pdx = plot['pdx'][plot['px'] - 0.5 * plot['pdx'] < 0]
        add2_x_pdy = plot['pdy'][plot['px'] - 0.5 * plot['pdx'] < 0]
        add2_x_field = plot[field][plot['px'] - 0.5 * plot['pdx'] < 0]
        add2_x_weight_field = plot['weight_field'][plot['px'] - 0.5 * plot['pdx'] < 0]

        # Cells hanging off the top edge.
        add_y_px = plot['px'][plot['py'] + 0.5 * plot['pdy'] > width[1]]
        add_y_py = plot['py'][plot['py'] + 0.5 * plot['pdy'] > width[1]]
        add_y_py -= width[1]
        add_y_pdx = plot['pdx'][plot['py'] + 0.5 * plot['pdy'] > width[1]]
        add_y_pdy = plot['pdy'][plot['py'] + 0.5 * plot['pdy'] > width[1]]
        add_y_field = plot[field][plot['py'] + 0.5 * plot['pdy'] > width[1]]
        add_y_weight_field = plot['weight_field'][plot['py'] + 0.5 * plot['pdy'] > width[1]]

        # Cells hanging off the bottom edge.
        add2_y_px = plot['px'][plot['py'] - 0.5 * plot['pdy'] < 0]
        add2_y_py = plot['py'][plot['py'] - 0.5 * plot['pdy'] < 0]
        add2_y_py += width[1]
        add2_y_pdx = plot['pdx'][plot['py'] - 0.5 * plot['pdy'] < 0]
        add2_y_pdy = plot['pdy'][plot['py'] - 0.5 * plot['pdy'] < 0]
        add2_y_field = plot[field][plot['py'] - 0.5 * plot['pdy'] < 0]
        add2_y_weight_field = plot['weight_field'][plot['py'] - 0.5 * plot['pdy'] < 0]

        # Add the hanging cells back to the projection data.
        plot.data['px'] = na.concatenate([plot['px'],add_x_px,add_y_px,add2_x_px,add2_y_px])
        plot.data['py'] = na.concatenate([plot['py'],add_x_py,add_y_py,add2_x_py,add2_y_py])
        plot.data['pdx'] = na.concatenate([plot['pdx'],add_x_pdx,add_y_pdx,add2_x_pdx,add2_y_pdx])
        plot.data['pdy'] = na.concatenate([plot['pdy'],add_x_pdy,add_y_pdy,add2_x_pdy,add2_y_pdy])
        plot.data[field] = na.concatenate([plot[field],add_x_field,add_y_field,add2_x_field,add2_y_field])
        plot.data['weight_field'] = na.concatenate([plot['weight_field'],
                                                    add_x_weight_field,add_y_weight_field,add2_x_weight_field,add2_y_weight_field])

        # Delete original copies of hanging cells.
        del add_x_px,add_y_px,add2_x_px,add2_y_px
        del add_x_py,add_y_py,add2_x_py,add2_y_py
        del add_x_pdx,add_y_pdx,add2_x_pdx,add2_y_pdx
        del add_x_pdy,add_y_pdy,add2_x_pdy,add2_y_pdy
        del add_x_field,add_y_field,add2_x_field,add2_y_field
        del add_x_weight_field,add_y_weight_field,add2_x_weight_field,add2_y_weight_field
