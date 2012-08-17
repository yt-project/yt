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
    particle_data,nstars = prepare_star_particles(pf,star_particle_type,fle=fle,fre=fre,
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

    return fle,fre,ile,ire,dd,nleaf,nstars

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

def prepare_octree(pf,ile,start_level=0,debug=True,dd=None,center=None):
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
    print 'first few entries :',refined[:12]
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
    if l<2:
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
    levels[pos.output_pos]  = level
    is_leaf = (child_grid_index==-1) and (level>0)
    refined[pos.refined_pos] = not is_leaf #True is oct, False is leaf
    pos.refined_pos+= 1 
    if is_leaf: #never subdivide if we are on a superlevel
        #then we have hit a leaf cell; write it out
        for field_index in range(grid.fields.shape[0]):
            output[pos.output_pos,field_index] = \
                    grid.fields[field_index,cell_index[0],cell_index[1],cell_index[2]]
        pos.output_pos+= 1 
    else:
        assert child_grid_index>-1
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
                child_ile = cell_index + na.array(vertex)*2**(-level)
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
    idxst = dd["particle_type"] == star_type

    #make sure we select more than a single particle
    assert na.sum(idxst)>0
    if pos is None:
        pos = na.array([dd["particle_position_%s" % ax]
                        for ax in 'xyz']).transpose()
    idx = idxst & na.all(pos>fle,axis=1) & na.all(pos<fre,axis=1)
    assert na.sum(idx)>0
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
    
    #make sure we have nonzero particle number
    assert pd_table.data.shape[0]>0
    return pd_table,na.sum(idx)


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

def generate_sunrise_cameraset_positions(pf,sim_center,cameraset=None,**kwargs):
    if cameraset is None:
        cameraset =cameraset_vertex 
    campos =[]
    names = []
    dd = pf.h.all_data()
    for name, (scene_pos,scene_up, scene_rot)  in cameraset.iteritems():
        kwargs['scene_position']=scene_pos
        kwargs['scene_up']=scene_up
        kwargs['scene_rot']=scene_rot
        kwargs['dd']=dd
        line = generate_sunrise_camera_position(pf,sim_center,**kwargs)
        campos += line,
        names += name,
    return names,campos     

def generate_sunrise_camera_position(pf,sim_center,sim_axis_short=None,sim_axis_long=None,
                                     sim_sphere_radius=None,sim_halo_radius=None,
                                     scene_position=[0.0,0.0,1.0],scene_distance=None,
                                     scene_up=[0.,0.,1.],scene_fov=None,scene_rot=True,
                                     dd=None):
    """Translate the simulation to center on sim_center, 
    then rotate such that sim_up is along the +z direction. Then we are in the 
    'scene' basis coordinates from which scene_up and scene_offset are defined.
    Then a position vector, direction vector, up vector and angular field of view
    are returned. The 3-vectors are in absolute physical kpc, not relative to the center.
    The angular field of view is in radians. The 10 numbers should match the inputs to
    camera_positions in Sunrise.
    """

    sim_center = na.array(sim_center)
    if sim_sphere_radius is None:
        sim_sphere_radius = 10.0/pf['kpc']
    if sim_axis_short is None:
        if dd is None:
            dd = pf.h.all_data()
        pos = na.array([dd["particle_position_%s"%i] for i in "xyz"]).T
        idx = na.sqrt(na.sum((pos-sim_center)**2.0,axis=1))<sim_sphere_radius
        mas = dd["particle_mass"]
        pos = pos[idx]
        mas = mas[idx]
        mo_inertia = position_moment(pos,mas)
        eigva, eigvc = linalg.eig(mo_inertia)
        #order into short, long axes
        order = eigva.real.argsort()
        ax_short,ax_med,ax_long = [ eigvc[:,order[i]] for i in (0,1,2)]
    else:
        ax_short = sim_axis_short
        ax_long  = sim_axis_long
    if sim_halo_radius is None:
        sim_halo_radius = 200.0/pf['kpc']
    if scene_distance is  None:
        scene_distance = 1e4/pf['kpc'] #this is how far the camera is from the target
    if scene_fov is None:
        radii = na.sqrt(na.sum((pos-sim_center)**2.0,axis=1))
        #idx= radii < sim_halo_radius*0.10
        #radii = radii[idx]
        #mass  = mas[idx] #copying mass into mas
        si = na.argsort(radii)
        radii = radii[si]
        mass  = mas[si]
        idx, = na.where(na.cumsum(mass)>mass.sum()/2.0)
        re = radii[idx[0]]
        scene_fov = 5*re
        scene_fov = max(scene_fov,3.0/pf['kpc']) #min size is 3kpc
        scene_fov = min(scene_fov,20.0/pf['kpc']) #max size is 3kpc
    #find rotation matrix
    angles=find_half_euler_angles(ax_short,ax_long)
    rotation  = euler_matrix(*angles)
    irotation = numpy.linalg.inv(rotation)
    axs = (ax_short,ax_med,ax_long)
    ax_rs,ax_rm,ax_rl = (matmul(rotation,ax) for ax in axs)
    axs = ([1,0,0],[0,1,0],[0,0,1])
    ax_is,ax_im,ax_il = (matmul(irotation,ax) for ax in axs)
    
    #rotate the camera
    if scene_rot :
        irotation = na.eye(3)
    sunrise_pos = matmul(irotation,na.array(scene_position)*scene_distance) #do NOT include sim center
    sunrise_up  = matmul(irotation,scene_up)
    sunrise_direction = -sunrise_pos
    sunrise_afov = 2.0*na.arctan((scene_fov/2.0)/scene_distance)#convert from distance FOV to angular

    #change to physical kpc
    sunrise_pos *= pf['kpc']
    sunrise_direction *= pf['kpc']
    return sunrise_pos,sunrise_direction,sunrise_up,sunrise_afov,scene_fov

def matmul(m, v):
    """Multiply a matrix times a set of vectors, or a single vector.
    My nPart x nDim convention leads to two transpositions, which is
    why this is hidden away in a function.  Note that if you try to
    use this to muliply two matricies, it will think that you're
    trying to multiply by a set of vectors and all hell will break
    loose."""    
    assert type(v) is not na.matrix
    v = na.asarray(v)
    m, vs = [na.asmatrix(a) for a in (m, v)]

    result = na.asarray(na.transpose(m * na.transpose(vs)))    
    if len(v.shape) == 1:
        return result[0]
    return result


def mag(vs):
    """Compute the norms of a set of vectors or a single vector."""
    vs = na.asarray(vs)
    if len(vs.shape) == 1:
        return na.sqrt( (vs**2).sum() )
    return na.sqrt( (vs**2).sum(axis=1) )

def mag2(vs):
    """Compute the norms of a set of vectors or a single vector."""
    vs = na.asarray(vs)
    if len(vs.shape) == 1:
        return (vs**2).sum()
    return (vs**2).sum(axis=1)


def position_moment(rs, ms=None, axes=None):
    """Find second position moment tensor.
    If axes is specified, weight by the elliptical radius (Allgood 2005)"""
    rs = na.asarray(rs)
    Npart, N = rs.shape
    if ms is None: ms = na.ones(Npart)
    else: ms = na.asarray(ms)    
    if axes is not None:
        axes = na.asarray(axes,dtype=float64)
        axes = axes/axes.max()
        norms2 = mag2(rs/axes)
    else:
        norms2 = na.ones(Npart)
    M = ms.sum()
    result = na.zeros((N,N))
    # matrix is symmetric, so only compute half of it then fill in the
    # other half
    for i in range(N):
        for j in range(i+1):
            result[i,j] = ( rs[:,i] * rs[:,j] * ms / norms2).sum() / M
        
    result = result + result.transpose() - na.identity(N)*result
    return result
    


def find_half_euler_angles(v,w,check=True):
    """Find the passive euler angles that will make v lie along the z
    axis and w lie along the x axis.  v and w are uncertain up to
    inversions (ie, eigenvectors) so this routine removes degeneracies
    associated with that

    (old) Calculate angles to bring a body into alignment with the
    coordinate system.  If v1 is the SHORTEST axis and v2 is the
    LONGEST axis, then this will return the angle (Euler angles) to
    make the long axis line up with the x axis and the short axis line
    up with the x (z) axis for the 2 (3) dimensional case."""
    # Make sure the vectors are normalized and orthogonal
    mag = lambda x: na.sqrt(na.sum(x**2.0))
    v = v/mag(v)
    w = w/mag(w)    
    if check:
        if abs((v*w).sum()) / (mag(v)*mag(w)) > 1e-5: raise ValueError

    # Break eigenvector scaling degeneracy by forcing it to have a positive
    # z component
    if v[2] < 0: v = -v
    phi,theta = find_euler_phi_theta(v)

    # Rotate w according to phi,theta and then break inversion
    # degeneracy by requiring that resulting vector has positive
    # x component
    w_prime = euler_passive(w,phi,theta,0.)
    if w_prime[0] < 0: w_prime = -w_prime
    # Now last Euler angle should just be this:
    psi = na.arctan2(w_prime[1],w_prime[0])
    return phi, theta, psi

def find_euler_phi_theta(v):
    """Find (passive) euler angles that will make v point in the z
    direction"""
    # Make sure the vector is normalized
    v = v/mag(v)
    theta = na.arccos(v[2])
    phi = na.arctan2(v[0],-v[1])
    return phi,theta

def euler_matrix(phi, the, psi):
    """Make an Euler transformation matrix"""
    cpsi=na.cos(psi)
    spsi=na.sin(psi)
    cphi=na.cos(phi)
    sphi=na.sin(phi)
    cthe=na.cos(the)
    sthe=na.sin(the)
    m = na.mat(na.zeros((3,3)))
    m[0,0] = cpsi*cphi - cthe*sphi*spsi
    m[0,1] = cpsi*sphi + cthe*cphi*spsi
    m[0,2] = spsi*sthe
    m[1,0] = -spsi*cphi - cthe*sphi*cpsi
    m[1,1] = -spsi*sphi + cthe*cphi*cpsi 
    m[1,2] = cpsi*sthe
    m[2,0] = sthe*sphi
    m[2,1] = -sthe*cphi
    m[2,2] = cthe
    return m

def euler_passive(v, phi, the, psi):
    """Passive Euler transform"""
    m = euler_matrix(phi, the, psi)
    return matmul(m,v)


#the format for these camerasets is name,up vector,camera location, 
#rotate to the galaxy's up direction?
cameraset_compass = collections.OrderedDict([
    ['top',([0.,0.,1.],[0.,-1.,0],True)], #up is north=+y
    ['bottom',([0.,0.,-1.],[0.,-1.,0.],True)],#up is north=+y
    ['north',([0.,1.,0.],[0.,0.,-1.],True)],#up is along z
    ['south',([0.,-1.,0.],[0.,0.,-1.],True)],#up is along z
    ['east',([1.,0.,0.],[0.,0.,-1.],True)],#up is along z
    ['west',([-1.,0.,0.],[0.,0.,-1.],True)],#up is along z
    ['top-north',([0.,0.7071,0.7071],[0., 0., -1.],True)],
    ['top-south',([0.,-0.7071,0.7071],[0., 0., -1.],True)],
    ['top-east',([ 0.7071,0.,0.7071],[0., 0., -1.],True)],
    ['top-west',([-0.7071,0.,0.7071],[0., 0., -1.],True)]
    ])

cameraset_vertex = collections.OrderedDict([
    ['top',([0.,0.,1.],[0.,-1.,0],True)], #up is north=+y
    ['north',([0.,1.,0.],[0.,0.,-1.],True)],#up is along z
    ['top-north',([0.,0.7071,0.7071],[0., 0., -1.],True)],
    ['Z',([0.,0.,1.],[0.,-1.,0],False)], #up is north=+y
    ['Y',([0.,1.,0.],[0.,0.,-1.],False)],#up is along z
    ['ZY',([0.,0.7071,0.7071],[0., 0., -1.],False)]
    ])

#up is 45deg down from z, towards north
#'bottom-north':([0.,0.7071,-0.7071],[0., 0., -1.])
#up is -45deg down from z, towards north

cameraset_ring = collections.OrderedDict()

segments = 20
for angle in na.linspace(0,360,segments):
    pos = [na.cos(angle),0.,na.sin(angle)]
    vc  = [na.cos(90-angle),0.,na.sin(90-angle)] 
    cameraset_ring['02i'%angle]=(pos,vc)
            


