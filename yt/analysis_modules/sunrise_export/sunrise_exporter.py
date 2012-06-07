"""
Code to export from yt to Sunrise

Author: Chris Moody <juxtaposicion@gmail.com>
Affiliation: UCSC
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

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

try:
    import pyfits
except ImportError: 
    pass

import time
import numpy as na

from yt.funcs import *
import yt.utilities.lib as amr_utils
from yt.data_objects.universal_fields import add_field
from yt.mods import *

debug = True

def export_to_sunrise(pf, fn, star_particle_type, dle, dre,**kwargs):
    r"""Convert the contents of a dataset to a FITS file format that Sunrise
    understands.

    This function will accept a parameter file, and from that parameter file
    construct a depth-first octree containing all of the data in the parameter
    file.  This octree will be written to a FITS file.  It will probably be
    quite big, so use this function with caution!  Sunrise is a tool for
    generating synthetic spectra, available at
    http://sunrise.googlecode.com/ .

    Parameters
    ----------
    pf : `StaticOutput`
        The parameter file to convert.
    fn : string
        The filename of the output FITS file.
    dle : The domain left edge to extract
    dre : The domain rght edge to extract
        Array format is (nx,ny,nz) where each element is floating point
        in unitary position units where 0 is leftmost edge and 1
        the rightmost. 
        

    Notes
    -----
    Note that the process of generating simulated images from Sunrise will
    require substantial user input; see the Sunrise wiki at
    http://sunrise.googlecode.com/ for more information.

    """
    
    #we must round the dle,dre to the nearest root grid cells
    ile,ire,super_level= round_nearest_edge(pf,dle,dre)
    super_level -= 1 #we're off by one (so we don't need a correction if we span 2 cells)
    fle,fre = ile*1.0/pf.domain_dimensions, ire*1.0/pf.domain_dimensions
    mylog.info("rounding specified region:")
    mylog.info("from [%1.5f %1.5f %1.5f]-[%1.5f %1.5f %1.5f]"%(tuple(dle)+tuple(dre)))
    mylog.info("to   [%07i %07i %07i]-[%07i %07i %07i]"%(tuple(ile)+tuple(ire)))
    mylog.info("to   [%1.5f %1.5f %1.5f]-[%1.5f %1.5f %1.5f]"%(tuple(fle)+tuple(fre)))


    #Create the refinement hilbert octree in GRIDSTRUCTURE
    #For every leaf (not-refined) cell we have a column n GRIDDATA
    #Include mass_gas, mass_metals, gas_temp_m, gas_teff_m, cell_volume, SFR
    #since the octree always starts with one cell, an our 0-level mesh
    #may have many cells, we must #create the octree region sitting 
    #ontop of the first mesh by providing a negative level
    output, refinement = prepare_octree(pf,ile,start_level=-super_level)

    #Create a list of the star particle properties in PARTICLE_DATA
    #Include ID, parent-ID, position, velocity, creation_mass, 
    #formation_time, mass, age_m, age_l, metallicity, L_bol
    particle_data = prepare_star_particles(pf,star_particle_type,fle=fle,fre=fre,**kwargs)

    create_fits_file(pf,fn, refinement,output,particle_data,fre,fle)

def prepare_octree(pf,ile,start_level=0):
    add_fields() #add the metal mass field that sunrise wants
    fields = ["CellMassMsun","TemperatureTimesCellMassMsun", 
              "MetalMass","CellVolumeCode"]
    
    #gather the field data from octs
    pbar = get_pbar("Retrieving field data",len(fields))
    field_data = [] 
    dd = pf.h.all_data()
    for fi,f in enumerate(fields):
        field_data += dd[f],
        pbar.update(fi)
    pbar.finish()
    del field_data

    #first we cast every cell as an oct
    #ngrids = na.max([g.id for g in pf._grids])
    grids = {}
    levels_all = {} 
    levels_finest = {}
    for l in range(100): 
        levels_finest[l]=0
        levels_all[l]=0
    pbar = get_pbar("Initializing octs ",len(pf.h.grids))
    for gi,g in enumerate(pf.h.grids):
        ff = na.array([g[f] for f in fields])
        og = amr_utils.OctreeGrid(
                g.child_index_mask.astype('int32'),
                ff.astype("float64"),
                g.LeftEdge.astype("float64"),
                g.ActiveDimensions.astype("int32"),
                na.ones(1,dtype="float64")*g.dds[0],
                g.Level,
                g.id)
        grids[g.id] = og
        #how many refinement cells will we have?
        #measure the 'volume' of each mesh, but many
        #cells do not exist. an overstimate
        levels_all[g.Level] += g.ActiveDimensions.prod()
        #how many leaves do we have?
        #this overestimates. a child of -1 means no child,
        #but that cell may still be expanded on a submesh because
        #(at least in ART) the meshes are inefficient.
        g.clear_data()
        pbar.update(gi)
    pbar.finish()
    
    #create the octree grid list
    oct_list =  amr_utils.OctreeGridList(grids)
    
    #initialize arrays to be passed to the recursion algo
    o_length = na.sum(levels_all.values())
    r_length = na.sum(levels_all.values())
    output   = na.zeros((o_length,len(fields)), dtype='float64')
    refined  = na.zeros(r_length, dtype='int32')
    levels   = na.zeros(r_length, dtype='int32')
    pos = position()
    hs       = hilbert_state()
    refined[0] = 1 #introduce the first cell as divided
    levels[0]  = start_level-1 #introduce the first cell as divided
    pos.refined_pos += 1
    RecurseOctreeDepthFirstHilbert(
            ile[0],ile[1],ile[2],
            pos,0, hs, 
            output,refined,levels,
            grids,
            start_level,
            #physical_center = (ile)*1.0/pf.domain_dimensions*pf['kpc'],
            physical_center = ile,
            #physical_width  = pf['kpc'])
            physical_width  = pf.domain_dimensions)
    #by time we get it here the 'current' position is actually 
    #for the next spot, so we're off by 1
    print 'refinement tree # of cells %i, # of leaves %i'%(pos.refined_pos,pos.output_pos) 
    output  = output[:pos.output_pos]
    refined = refined[:pos.refined_pos] 
    levels = levels[:pos.refined_pos] 
    return output,refined

def print_row(level,ple,pre,pw,pc,hs):
    print level, 
    print '%1.5f %1.5f %1.5f '%tuple(ple*pw-pc),
    print '%1.5f %1.5f %1.5f '%tuple(pre*pw-pc),
    print hs.dim, hs.sgn

def print_child(level,grid,i,j,k,pw,pc):
    ple = (grid.left_edges+na.array([i,j,k])*grid.dx)*pw-pc #parent LE 
    pre = (grid.left_edges+na.array([i+1,j+1,k+1])*grid.dx)*pw-pc #parent RE 
    print level, 
    print '%1.5f %1.5f %1.5f '%tuple(ple),
    print '%1.5f %1.5f %1.5f '%tuple(pre)

def RecurseOctreeDepthFirstHilbert(xi,yi,zi,
                            curpos, gi, 
                            hs,
                            output,
                            refined,
                            levels,
                            grids,
                            level,
                            physical_center=None,
                            physical_width=None,
                            printr=False):
    grid = grids[gi]
    m = 2**(-level-1) if level < 0 else 1
    ple = grid.left_edges+na.array([xi,yi,zi])*grid.dx #parent LE
    pre = ple+grid.dx*m
    if printr:
        print_row(level,ple,pre,physical_width,physical_center,hs)

    #here we go over the 8 octants
    #in general however, a mesh cell on this level
    #may have more than 8 children on the next level
    #so we find the int float center (cxyz) of each child cell
    # and from that find the child cell indices
    for iv, (vertex,hs_child) in enumerate(hs):
        #print ' '*(level+3), level,iv, vertex,curpos.refined_pos,curpos.output_pos,
        #negative level indicates that we need to build a super-octree
        if level < 0: 
            #print ' '
            #we are not on the root grid yet, but this is 
            #how many equivalent root grid cells we would have
            #level -1 means our oct grid's children are the same size
            #as the root grid (hence the -level-1)
            dx = 2**(-level-1) #this is the child width 
            i,j,k = xi+vertex[0]*dx,yi+vertex[1]*dx,zi+vertex[2]*dx
            #we always refine the negative levels
            refined[curpos.refined_pos] = 1
            levels[curpos.refined_pos] = level
            curpos.refined_pos += 1
            RecurseOctreeDepthFirstHilbert(i, j, k,
                                curpos, 0, hs_child, output, refined, levels, grids,
                                level+1,
                                physical_center=physical_center,
                                physical_width=physical_width,)
        else:
            i,j,k = xi+vertex[0],yi+vertex[1],zi+vertex[2]
            ci = grid.child_indices[i,j,k] #is this oct subdivided?
            if ci == -1:
                for fi in range(grid.fields.shape[0]):
                    output[curpos.output_pos,fi] = grid.fields[fi,i,j,k]
                refined[curpos.refined_pos] = 0
                levels[curpos.refined_pos] = level
                curpos.output_pos += 1 #position updated after write
                curpos.refined_pos += 1
                if printr:
                    print_child(level+1,grid,i,j,k,physical_width,physical_center)
            else:
                cx = (grid.left_edges[0] + i*grid.dx[0]) #floating le of the child
                cy = (grid.left_edges[1] + j*grid.dx[0])
                cz = (grid.left_edges[2] + k*grid.dx[0])
                refined[curpos.refined_pos] = 1
                levels[curpos.refined_pos] = level
                curpos.refined_pos += 1 #position updated after write
                child_grid = grids[ci]
                child_dx = child_grid.dx[0]
                child_leftedges = child_grid.left_edges
                child_i = int((cx - child_leftedges[0])/child_dx)
                child_j = int((cy - child_leftedges[1])/child_dx)
                child_k = int((cz - child_leftedges[2])/child_dx)
                RecurseOctreeDepthFirstHilbert(child_i, child_j, child_k,
                                    curpos, ci, hs_child, output, refined, levels, grids,
                                    level+1,
                                    physical_center=physical_center,
                                    physical_width=physical_width)

def create_fits_file(pf,fn, refined,output,particle_data,fre,fle):

    #first create the grid structure
    structure = pyfits.Column("structure", format="B", array=refined.astype("bool"))
    cols = pyfits.ColDefs([structure])
    st_table = pyfits.new_table(cols)
    st_table.name = "GRIDSTRUCTURE"
    st_table.header.update("hierarch lengthunit", "kpc", comment="Length unit for grid")
    fdx = fre-fle
    for i,a in enumerate('xyz'):
        st_table.header.update("min%s" % a, fle[i] * pf['kpc'])
        st_table.header.update("max%s" % a, fre[i] * pf['kpc'])
        #st_table.header.update("min%s" % a, 0) #WARNING: this is for debugging
        #st_table.header.update("max%s" % a, 2) #
        st_table.header.update("n%s" % a, fdx[i])
        st_table.header.update("subdiv%s" % a, 2)
    st_table.header.update("subdivtp", "OCTREE", "Type of grid subdivision")

    #not the hydro grid data
    fields = ["CellMassMsun","TemperatureTimesCellMassMsun", 
              "MetalMass","CellVolumeCode"]
    fd = {}
    for i,f in enumerate(fields): 
        fd[f]=output[:,i]
    del output
    col_list = []
    size = fd["CellMassMsun"].size
    tm = fd["CellMassMsun"].sum()
    col_list.append(pyfits.Column("mass_gas", format='D',
                    array=fd['CellMassMsun'], unit="Msun"))
    col_list.append(pyfits.Column("mass_metals", format='D',
                    array=fd['MetalMass'], unit="Msun"))
    # col_list.append(pyfits.Column("mass_stars", format='D',
    #                 array=na.zeros(size,dtype='D'),unit="Msun"))
    # col_list.append(pyfits.Column("mass_stellar_metals", format='D',
    #                 array=na.zeros(size,dtype='D'),unit="Msun"))
    # col_list.append(pyfits.Column("age_m", format='D',
    #                 array=na.zeros(size,dtype='D'),unit="yr*Msun"))
    # col_list.append(pyfits.Column("age_l", format='D',
    #                 array=na.zeros(size,dtype='D'),unit="yr*Msun"))
    # col_list.append(pyfits.Column("L_bol", format='D',
    #                 array=na.zeros(size,dtype='D')))
    # col_list.append(pyfits.Column("L_lambda", format='D',
    #                 array=na.zeros(size,dtype='D')))
    # The units for gas_temp are really K*Msun. For older Sunrise versions
    # you must set the unit to just K  
    col_list.append(pyfits.Column("gas_temp_m", format='D',
                    array=fd['TemperatureTimesCellMassMsun'], unit="K*Msun"))
    col_list.append(pyfits.Column("gas_teff_m", format='D',
                    array=fd['TemperatureTimesCellMassMsun'], unit="K*Msun"))
    col_list.append(pyfits.Column("cell_volume", format='D',
                    array=fd['CellVolumeCode'].astype('float64')*pf['kpc']**3.0,
                    unit="kpc^3"))
    col_list.append(pyfits.Column("SFR", format='D',
                    array=na.zeros(size, dtype='D')))
    cols = pyfits.ColDefs(col_list)
    mg_table = pyfits.new_table(cols)
    mg_table.header.update("M_g_tot", tm)
    mg_table.header.update("timeunit", "yr")
    mg_table.header.update("tempunit", "K")
    mg_table.name = "GRIDDATA"

    # Add a dummy Primary; might be a better way to do this!
    col_list = [pyfits.Column("dummy", format="F", array=na.zeros(1, dtype='float32'))]
    cols = pyfits.ColDefs(col_list)
    md_table = pyfits.new_table(cols)
    md_table.header.update("snaptime", pf.current_time*pf['years'])
    md_table.name = "YT"

    phdu = pyfits.PrimaryHDU()
    phdu.header.update('nbodycod','yt')
    hls = [phdu, st_table, mg_table,md_table]
    hls.append(particle_data)
    hdus = pyfits.HDUList(hls)
    hdus.writeto(fn, clobber=True)

def nearest_power(x):
    #round to the nearest power of 2
    x-=1
    x |= x >> 1
    x |= x >> 2 
    x |= x >> 4
    x |= x >> 8
    x |= x >> 16
    x+=1 
    return x

def round_nearest_edge(pf,dle,dre):
    dds = pf.domain_dimensions
    ile = na.floor(dle*dds).astype('int')
    ire = na.ceil(dre*dds).astype('int') 
    
    #this is the number of cells the super octree needs to expand to
    #must round to the nearest power of 2
    width = na.max(ire-ile)
    width = nearest_power(width)
    
    maxlevel = na.rint(na.log2(width)).astype('int')
    return ile,ire,maxlevel

def prepare_star_particles(pf,star_type,pos=None,vel=None, age=None,
                          creation_time=None,initial_mass=None,
                          current_mass=None,metallicity=None,
                          radius = None,
                          fle=[0.,0.,0.],fre=[1.,1.,1.]):
    dd = pf.h.all_data()
    idx = dd["particle_type"] == star_type
    if pos is None:
        pos = na.array([dd["particle_position_%s" % ax]
                        for ax in 'xyz']).transpose()
    idx = idx & na.all(pos>fle,axis=1) & na.all(pos<fre,axis=1)
    pos = pos[idx]*pf['kpc'] #unitary units -> kpc
    if age is None:
        age = dd["particle_age"][idx]*pf['years'] # seconds->years
    if vel is None:
        vel = na.array([dd["particle_velocity_%s" % ax][idx]
                        for ax in 'xyz']).transpose()
        # Velocity is cm/s, we want it to be kpc/yr
        #vel *= (pf["kpc"]/pf["cm"]) / (365*24*3600.)
        vel *= 1.02268944e-14 
    if initial_mass is None:
        #in solar masses
        initial_mass = dd["particle_mass_initial"][idx]*pf['Msun']
    if current_mass is None:
        #in solar masses
        current_mass = dd["particle_mass"][idx]*pf['Msun']
    if metallicity is None:
        #this should be in dimensionless units, metals mass / particle mass
        metallicity = dd["particle_metallicity"][idx]
    if radius is None:
        radius = initial_mass*0.0+10.0/1000.0 #10pc radius

    formation_time = pf.current_time-age
    #create every column
    col_list = []
    col_list.append(pyfits.Column("ID", format="I", array=na.arange(current_mass.size)))
    col_list.append(pyfits.Column("parent_ID", format="I", array=na.arange(current_mass.size)))
    col_list.append(pyfits.Column("position", format="3D", array=pos, unit="kpc"))
    col_list.append(pyfits.Column("velocity", format="3D", array=vel, unit="kpc/yr"))
    col_list.append(pyfits.Column("creation_mass", format="D", array=initial_mass, unit="Msun"))
    col_list.append(pyfits.Column("formation_time", format="D", array=formation_time, unit="yr"))
    col_list.append(pyfits.Column("radius", format="D", array=radius, unit="kpc"))
    col_list.append(pyfits.Column("mass", format="D", array=current_mass, unit="Msun"))
    col_list.append(pyfits.Column("age_m", format="D", array=age))
    col_list.append(pyfits.Column("age_l", format="D", array=age))
    #For particles, Sunrise takes 
    #the dimensionless metallicity, not the mass of the metals
    col_list.append(pyfits.Column("metallicity", format="D",
        array=metallicity,unit="Msun")) 
    col_list.append(pyfits.Column("L_bol", format="D",
        array=na.zeros(current_mass.size)))
    
    #make the table
    cols = pyfits.ColDefs(col_list)
    pd_table = pyfits.new_table(cols)
    pd_table.name = "PARTICLEDATA"
    return pd_table


def add_fields():
    """Add three Eulerian fields Sunrise uses"""
    def _MetalMass(field, data):
        return data["Metal_Density"] * data["CellVolume"]
        
    def _convMetalMass(data):
        return 1.0/1.989e33
    
    add_field("MetalMass", function=_MetalMass,
              convert_function=_convMetalMass)

    def _initial_mass_cen_ostriker(field, data):
        # SFR in a cell. This assumes stars were created by the Cen & Ostriker algorithm
        # Check Grid_AddToDiskProfile.C and star_maker7.src
        star_mass_ejection_fraction = data.pf.get_parameter("StarMassEjectionFraction",float)
        star_maker_minimum_dynamical_time = 3e6 # years, which will get divided out
        dtForSFR = star_maker_minimum_dynamical_time / data.pf["years"]
        xv1 = ((data.pf["InitialTime"] - data["creation_time"])
                / data["dynamical_time"])
        xv2 = ((data.pf["InitialTime"] + dtForSFR - data["creation_time"])
                / data["dynamical_time"])
        denom = (1.0 - star_mass_ejection_fraction * (1.0 - (1.0 + xv1)*na.exp(-xv1)))
        minitial = data["ParticleMassMsun"] / denom
        return minitial

    add_field("InitialMassCenOstriker", function=_initial_mass_cen_ostriker)

    def _temp_times_mass(field, data):
        return data["Temperature"]*data["CellMassMsun"]
    add_field("TemperatureTimesCellMassMsun", function=_temp_times_mass)

class position:
    def __init__(self):
        self.output_pos = 0
        self.refined_pos = 0

class hilbert_state():
    def __init__(self,dim=None,sgn=None,octant=None):
        if dim is None: dim = [0,1,2]
        if sgn is None: sgn = [1,1,1]
        if octant is None: octant = 5
        self.dim = dim
        self.sgn = sgn
        self.octant = octant
    def flip(self,i):
        self.sgn[i]*=-1
    def swap(self,i,j):
        temp = self.dim[i]
        self.dim[i]=self.dim[j]
        self.dim[j]=temp
        axis = self.sgn[i]
        self.sgn[i] = self.sgn[j]
        self.sgn[j] = axis
    def reorder(self,i,j,k):
        ndim = [self.dim[i],self.dim[j],self.dim[k]] 
        nsgn = [self.sgn[i],self.sgn[j],self.sgn[k]]
        self.dim = ndim
        self.sgn = nsgn
    def copy(self):
        return hilbert_state([self.dim[0],self.dim[1],self.dim[2]],
                             [self.sgn[0],self.sgn[1],self.sgn[2]],
                             self.octant)
    def descend(self,o):
        child = self.copy()
        child.octant = o
        if o==0:
            child.swap(0,2)
        elif o==1:
            child.swap(1,2)
        elif o==2:
            pass
        elif o==3:
            child.flip(0)
            child.flip(2)
            child.reorder(2,0,1)
        elif o==4:
            child.flip(0)
            child.flip(1)
            child.reorder(2,0,1)
        elif o==5:
            pass
        elif o==6:
            child.flip(1)
            child.flip(2)
            child.swap(1,2)
        elif o==7:
            child.flip(0)
            child.flip(2)
            child.swap(0,2)
        return child

    def __iter__(self):
        vertex = [0,0,0]
        j=0
        for i in range(3):
            vertex[self.dim[i]] = 0 if self.sgn[i]>0 else 1
        yield vertex, self.descend(j)
        vertex[self.dim[0]] += self.sgn[0]
        j+=1
        yield vertex, self.descend(j)
        vertex[self.dim[1]] += self.sgn[1] 
        j+=1
        yield vertex, self.descend(j)
        vertex[self.dim[0]] -= self.sgn[0] 
        j+=1
        yield vertex, self.descend(j)
        vertex[self.dim[2]] += self.sgn[2] 
        j+=1
        yield vertex, self.descend(j)
        vertex[self.dim[0]] += self.sgn[0] 
        j+=1
        yield vertex, self.descend(j)
        vertex[self.dim[1]] -= self.sgn[1] 
        j+=1
        yield vertex, self.descend(j)
        vertex[self.dim[0]] -= self.sgn[0] 
        j+=1
        yield vertex, self.descend(j)





