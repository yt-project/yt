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
import numpy.linalg as linalg
import collections

from yt.funcs import *
import yt.utilities.lib as amr_utils
from yt.data_objects.universal_fields import add_field
from yt.mods import *

debug = True

def export_to_sunrise(pf, fn, star_particle_type, fc, fwidth, ncells_wide=None,
        debug=False,dd=None,**kwargs):
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
    pf      : `StaticOutput`
                The parameter file to convert.
    fn      : string
                The filename of the output FITS file.
    fc      : array
                The center of the extraction region
    fwidth  : array  
                Ensure this radius around the center is enclosed
        Array format is (nx,ny,nz) where each element is floating point
        in unitary position units where 0 is leftmost edge and 1
        the rightmost. 
        

    Notes
    -----
    Note that the process of generating simulated images from Sunrise will
    require substantial user input; see the Sunrise wiki at
    http://sunrise.googlecode.com/ for more information.

    """

    fc = na.array(fc)
    fwidth = na.array(fwidth)
    
    #we must round the dle,dre to the nearest root grid cells
    ile,ire,super_level,ncells_wide= \
            round_ncells_wide(pf.domain_dimensions,fc-fwidth,fc+fwidth,nwide=ncells_wide)

    assert na.all((ile-ire)==(ile-ire)[0])
    mylog.info("rounding specified region:")
    mylog.info("from [%1.5f %1.5f %1.5f]-[%1.5f %1.5f %1.5f]"%(tuple(fc-fwidth)+tuple(fc+fwidth)))
    mylog.info("to   [%07i %07i %07i]-[%07i %07i %07i]"%(tuple(ile)+tuple(ire)))
    fle,fre = ile*1.0/pf.domain_dimensions, ire*1.0/pf.domain_dimensions
    mylog.info("to   [%1.5f %1.5f %1.5f]-[%1.5f %1.5f %1.5f]"%(tuple(fle)+tuple(fre)))

    #Create a list of the star particle properties in PARTICLE_DATA
    #Include ID, parent-ID, position, velocity, creation_mass, 
    #formation_time, mass, age_m, age_l, metallicity, L_bol
    particle_data = prepare_star_particles(pf,star_particle_type,fle=fle,fre=fre,
                                           dd=dd,**kwargs)

    #Create the refinement hilbert octree in GRIDSTRUCTURE
    #For every leaf (not-refined) cell we have a column n GRIDDATA
    #Include mass_gas, mass_metals, gas_temp_m, gas_teff_m, cell_volume, SFR
    #since the octree always starts with one cell, an our 0-level mesh
    #may have many cells, we must create the octree region sitting 
    #ontop of the first mesh by providing a negative level
    output, refinement,dd,nleaf = prepare_octree(pf,ile,start_level=super_level,
            debug=debug,dd=dd,center=fc)

    create_fits_file(pf,fn, refinement,output,particle_data,fle,fre)

    return fle,fre,ile,ire,dd,nleaf

def export_to_sunrise_from_halolist(pf,fni,star_particle_type,
                                        halo_list,domains_list=None,**kwargs):
    """
    Using the center of mass and the virial radius
    for a halo, calculate the regions to extract for sunrise.
    The regions are defined on the root grid, and so individual
    octs may span a large range encompassing many halos
    and subhalos. Instead of repeating the oct extraction for each
    halo, arrange halos such that we only calculate what we need to.

    Parameters
    ----------
    pf : `StaticOutput`
        The parameter file to convert. We use the root grid to specify the domain.
    fni : string
        The filename of the output FITS file, but depends on the domain. The
        dle and dre are appended to the name.
    particle_type : int
        The particle index for stars
    halo_list : list of halo objects
        The halo list objects must have halo.CoM and halo.Rvir,
        both of which are assumed to be in unitary length units.
    frvir (optional) : float
        Ensure that CoM +/- frvir*Rvir is contained within each domain
    domains_list (optiona): dict of halos
        Organize halos into a dict of domains. Keys are DLE/DRE tuple
        values are a list of halos
    """
    dn = pf.domain_dimensions
    if domains_list is None:
        domains_list = domains_from_halos(pf,halo_list,**kwargs)
    if fni.endswith('.fits'):
        fni = fni.replace('.fits','')

    ndomains_finished = 0
    for (num_halos, domain, halos) in domains_list:
        dle,dre = domain
        print 'exporting: '
        print "[%03i %03i %03i] -"%tuple(dle),
        print "[%03i %03i %03i] "%tuple(dre),
        print " with %i halos"%num_halos
        dle,dre = domain
        dle, dre = na.array(dle),na.array(dre)
        fn = fni 
        fn += "%03i_%03i_%03i-"%tuple(dle)
        fn += "%03i_%03i_%03i"%tuple(dre)
        fnf = fn + '.fits'
        fnt = fn + '.halos'
        if os.path.exists(fnt):
            os.remove(fnt)
        fh = open(fnt,'w')
        for halo in halos:
            fh.write("%i "%halo.ID)
            fh.write("%6.6e "%(halo.CoM[0]*pf['kpc']))
            fh.write("%6.6e "%(halo.CoM[1]*pf['kpc']))
            fh.write("%6.6e "%(halo.CoM[2]*pf['kpc']))
            fh.write("%6.6e "%(halo.Mvir))
            fh.write("%6.6e \n"%(halo.Rvir*pf['kpc']))
        fh.close()
        export_to_sunrise(pf, fnf, star_particle_type, dle*1.0/dn, dre*1.0/dn)
        ndomains_finished +=1

def domains_from_halos(pf,halo_list,frvir=0.15):
    domains = {}
    dn = pf.domain_dimensions
    for halo in halo_list:
        fle, fre = halo.CoM-frvir*halo.Rvir,halo.CoM+frvir*halo.Rvir
        dle,dre = na.floor(fle*dn), na.ceil(fre*dn)
        dle,dre = tuple(dle.astype('int')),tuple(dre.astype('int'))
        if (dle,dre) in domains.keys():
            domains[(dle,dre)] += halo,
        else:
            domains[(dle,dre)] = [halo,]
    #for niceness, let's process the domains in order of 
    #the one with the most halos
    domains_list = [(len(v),k,v) for k,v in domains.iteritems()]
    domains_list.sort() 
    domains_list.reverse() #we want the most populated domains first
    domains_limits = [d[1] for d in domains_list]
    domains_halos  = [d[2] for d in domains_list]
    return domains_list

def prepare_octree(pf,ile,start_level=0,debug=False,dd=None,center=None):
    add_fields() #add the metal mass field that sunrise wants
    fields = ["CellMassMsun","TemperatureTimesCellMassMsun", 
              "MetalMass","CellVolumeCode"]
    
    #gather the field data from octs
    pbar = get_pbar("Retrieving field data",len(fields))
    field_data = [] 
    if dd is None:
        #we keep passing dd around to not regenerate the data all the time
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
    #oct_list =  amr_utils.OctreeGridList(grids)
    
    #initialize arrays to be passed to the recursion algo
    o_length = na.sum(levels_all.values())
    r_length = na.sum(levels_all.values())
    output   = na.zeros((o_length,len(fields)), dtype='float64')
    refined  = na.zeros(r_length, dtype='int32')
    levels   = na.zeros(r_length, dtype='int32')
    pos = position()
    hs       = hilbert_state()
    start_time = time.time()
    if debug:
        if center is not None: 
            c = center*pf['kpc']
        else:
            c = ile*1.0/pf.domain_dimensions*pf['kpc']
        printing = lambda x: print_oct(x,pf['kpc'],c)
    else:
        printing = None
    pbar = get_pbar("Building Hilbert DFO octree",len(refined))
    RecurseOctreeDepthFirstHilbert(
            ile,
            pos,
            grids[0], #we always start on the root grid
            hs, 
            output,refined,levels,
            grids,
            start_level,
            debug=printing,
            tracker=pbar)
    pbar.finish()
    #by time we get it here the 'current' position is actually 
    #for the next spot, so we're off by 1
    print 'took %1.2e seconds'%(time.time()-start_time)
    print 'refinement tree # of cells %i, # of leaves %i'%(pos.refined_pos,pos.output_pos) 
    output  = output[:pos.output_pos]
    refined = refined[:pos.refined_pos] 
    levels = levels[:pos.refined_pos] 
    return output,refined,dd,pos.refined_pos

def print_oct(data,nd=None,nc=None):
    ci = data['cell_index']
    l  = data['level']
    g  = data['grid']
    fle = g.left_edges+g.dx*ci
    fre = g.left_edges+g.dx*(ci+1)
    if nd is not None:
        fle *= nd
        fre *= nd
        if nc is not None:
            fle -= nc
            fre -= nc
    txt  = '%1i '
    txt += '%1.3f '*3+'- '
    txt += '%1.3f '*3
    print txt%((l,)+tuple(fle)+tuple(fre))

def RecurseOctreeDepthFirstHilbert(cell_index, #integer (rep as a float) on the grids[grid_index]
                            pos, #the output hydro data position and refinement position
                            grid,  #grid that this oct lives on (not its children)
                            hilbert,  #the hilbert state
                            output, #holds the hydro data
                            refined, #holds the refinement status  of Octs, 0s and 1s
                            levels, #For a given Oct, what is the level
                            grids, #list of all patch grids available to us
                            level, #starting level of the oct (not the children)
                            debug=None,tracker=True):
    if tracker is not None:
        if pos.refined_pos%1000 == 500 : tracker.update(pos.refined_pos)
    if debug is not None: 
        debug(vars())
    child_grid_index = grid.child_indices[cell_index[0],cell_index[1],cell_index[2]]
    #record the refinement state
    refined[pos.refined_pos] = child_grid_index!=-1
    levels[pos.output_pos]  = level
    pos.refined_pos+= 1 
    if child_grid_index == -1 and level>=0: #never subdivide if we are on a superlevel
        #then we have hit a leaf cell; write it out
        for field_index in range(grid.fields.shape[0]):
            output[pos.output_pos,field_index] = \
                    grid.fields[field_index,cell_index[0],cell_index[1],cell_index[2]]
        pos.output_pos+= 1 
    else:
        #find the grid we descend into
        #then find the eight cells we break up into
        subgrid = grids[child_grid_index]
        #calculate the floating point LE of the children
        #then translate onto the subgrid integer index 
        parent_fle  = grid.left_edges + cell_index*grid.dx
        subgrid_ile = na.floor((parent_fle - subgrid.left_edges)/subgrid.dx)
        for i, (vertex,hilbert_child) in enumerate(hilbert):
            #vertex is a combination of three 0s and 1s to 
            #denote each of the 8 octs
            if level < 0:
                subgrid = grid #we don't actually descend if we're a superlevel
                child_ile = cell_index + vertex*2**(-level)
            else:
                child_ile = subgrid_ile+na.array(vertex)
                child_ile = child_ile.astype('int')
            RecurseOctreeDepthFirstHilbert(child_ile,pos,
                    subgrid,hilbert_child,output,refined,levels,grids,level+1,
                    debug=debug,tracker=tracker)



def create_fits_file(pf,fn, refined,output,particle_data,fle,fre):

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

def round_ncells_wide(dds,fle,fre,nwide=None):
    fc = (fle+fre)/2.0
    assert na.all(fle < fc)
    assert na.all(fre > fc)
    ic = na.rint(fc*dds) #nearest vertex to the center
    ile,ire = ic.astype('int'),ic.astype('int')
    cfle,cfre = fc.copy(),fc.copy()
    idx = na.array([0,0,0]) #just a random non-equal array
    width = 0.0
    if nwide is None:
        #expand until borders are included and
        #we have an equaly-sized, non-zero box
        idxq,out=False,True
        while not out or not idxq:
            cfle,cfre = fc-width, fc+width
            ile = na.rint(cfle*dds).astype('int')
            ire = na.rint(cfre*dds).astype('int')
            idx = ire-ile
            width += 0.1/dds
            #quit if idxq is true:
            idxq = idx[0]>0 and na.all(idx==idx[0])
            out  = na.all(fle>cfle) and na.all(fre<cfre) 
            assert width[0] < 1.1 #can't go larger than the simulation volume
        nwide = idx[0]
    else:
        #expand until we are nwide cells span
        while not na.all(idx==nwide):
            assert na.any(idx<=nwide)
            cfle,cfre = fc-width, fc+width
            ile = na.rint(cfle*dds).astype('int')
            ire = na.rint(cfre*dds).astype('int')
            idx = ire-ile
            width += 1e-2*1.0/dds
    assert na.all(idx==nwide)
    assert idx[0]>0
    maxlevel = -na.rint(na.log2(nwide)).astype('int')
    assert abs(na.log2(nwide)-na.rint(na.log2(nwide)))<1e-5 #nwide should be a power of 2
    return ile,ire,maxlevel,nwide

def round_nearest_edge(pf,fle,fre):
    dds = pf.domain_dimensions
    ile = na.floor(fle*dds).astype('int')
    ire = na.ceil(fre*dds).astype('int') 
    
    #this is the number of cells the super octree needs to expand to
    #must round to the nearest power of 2
    width = na.max(ire-ile)
    width = nearest_power(width)
    
    maxlevel = -na.rint(na.log2(width)).astype('int')
    return ile,ire,maxlevel

def prepare_star_particles(pf,star_type,pos=None,vel=None, age=None,
                          creation_time=None,initial_mass=None,
                          current_mass=None,metallicity=None,
                          radius = None,
                          fle=[0.,0.,0.],fre=[1.,1.,1.],
                          dd=None):
    if dd is None:
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
        #metallicity *=0.0198
        #print 'WARNING: multiplying metallicirt by 0.0198'
    if radius is None:
        radius = initial_mass*0.0+10.0/1000.0 #10pc radius
    formation_time = pf.current_time*pf['years']-age
    #create every column
    col_list = []
    col_list.append(pyfits.Column("ID", format="J", array=na.arange(current_mass.size).astype('int32')))
    col_list.append(pyfits.Column("parent_ID", format="J", array=na.arange(current_mass.size).astype('int32')))
    col_list.append(pyfits.Column("position", format="3D", array=pos, unit="kpc"))
    col_list.append(pyfits.Column("velocity", format="3D", array=vel, unit="kpc/yr"))
    col_list.append(pyfits.Column("creation_mass", format="D", array=initial_mass, unit="Msun"))
    col_list.append(pyfits.Column("formation_time", format="D", array=formation_time, unit="yr"))
    col_list.append(pyfits.Column("radius", format="D", array=radius, unit="kpc"))
    col_list.append(pyfits.Column("mass", format="D", array=current_mass, unit="Msun"))
    col_list.append(pyfits.Column("age", format="D", array=age,unit='yr'))
    #col_list.append(pyfits.Column("age_l", format="D", array=age, unit = 'yr'))
    #For particles, Sunrise takes 
    #the dimensionless metallicity, not the mass of the metals
    col_list.append(pyfits.Column("metallicity", format="D",
        array=metallicity,unit="Msun")) 
    #col_list.append(pyfits.Column("L_bol", format="D",
    #    array=na.zeros(current_mass.size)))
    
    #make the table
    cols = pyfits.ColDefs(col_list)
    pd_table = pyfits.new_table(cols)
    pd_table.name = "PARTICLEDATA"
    return pd_table


def add_fields():
    """Add three Eulerian fields Sunrise uses"""
    def _MetalMass(field, data):
        return data["Metallicity"] * data["CellMassMsun"]
        
    def _convMetalMass(data):
        return 1.0
    
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





