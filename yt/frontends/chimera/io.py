r"""
    Chimera-specific IO functions
    
    
    
    """

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import sys
import numpy as np
from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.on_demand_imports import _h5py as h5py
from yt.utilities.file_handler import \
    HDF5FileHandler
from yt.frontends.chimera.data_structures import _find_files




class ChimeraIOHandler(BaseIOHandler):
    _particle_reader = False
    _dataset_type = 'chimera'
    
    def __init__(self, ds):
        super(ChimeraIOHandler, self).__init__(ds)
        self._handle = ds._handle
        self.filename = ds.filename
    
    def _read_particle_coords(self, chunks, ptf):
        # This needs to *yield* a series of tuples of (ptype, (x, y, z)).
        # chunks is a list of chunks, and ptf is a dict where the keys are
        # ptypes and the values are lists of fields.
        pass
    
    def _read_particle_fields(self, chunks, ptf, selector):
        # This gets called after the arrays have been allocated.  It needs to
        # yield ((ptype, field), data) where data is the masked results of
        # reading ptype, field and applying the selector to the data read in.
        # Selector objects have a .select_points(x,y,z) that returns a mask, so
        # you need to do your masking here.
        pass
    
    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        nodal_fields = []
        for field in fields:
            finfo = self.ds.field_info[field]
            nodal_flag = finfo.nodal_flag
            if np.any(nodal_flag):
                num_nodes = 2**sum(nodal_flag)
                rv[field] = np.empty((size, num_nodes), dtype="=f8")
                nodal_fields.append(field)
            else:
                rv[field] = np.empty(size, dtype="=f8")
        ind = 0
        for field, mesh, data in self.io_iter(chunks, fields):
            if data is None:
                continue
            else:
                ind += mesh.select(selector, data.flatten(), rv[field], ind)  # caches
        return rv
    
    def io_iter(self, chunks, fields):
        for n, chunk in enumerate(chunks):
            file = _find_files(self.filename)
            if any('grid_2' in f for f in file):
                yy = True
            else:
                yy = False
            with h5py.File(file[n],"r") as f:
                #Generates mask according to the "ongrid_mask" variable
                m = int(file[n][-5:-3])-1
                k = f["fluid"]["entropy"].shape[0]
                mask_0 = f["mesh"]["ongrid_mask"][k*m:k*(m+1),:]
                #mask_0 = np.random.randint(0,2, size = (15,90)).astype('float')
                
                if (f["mesh"]["array_dimensions"][2] > 1):
                    nrd = f["mesh"]["array_dimensions"][0]-2
                else:
                    nrd = f["mesh"]["array_dimensions"][0]
            
                mask = np.repeat(mask_0[:, :, np.newaxis], nrd, axis=2).transpose()
                for field in fields: #Reads data by locating subheading
                    ftype, fname = field
                    specials = ('abar', 'e_rms_1','e_rms_2','e_rms_3','e_rms_4',
                                'lumin_1','lumin_2','lumin_3','lumin_4',
                                'num_lumin_1','num_lumin_2','num_lumin_3','num_lumin_4','shock')
                    a_name_2 = [i.decode("utf-8") for i in f["abundance"]["a_name"]]
                    if fname not in specials:
                        if fname in f["fluid"]:
                            ds = f["fluid"]["{}".format(fname)]
                        elif fname in f["abundance"]:
                            ds = f["abundance"]["{}".format(fname)]
                        elif fname in a_name_2:
                            ind_xn = a_name_2.index(fname)
                            ds = f["abundance"]["xn_c"][:,:,:,ind_xn]
                        else:
                            sys.exit('Error: Invalid field name')
                        dat_1 = ds[:,:,:].transpose()
                                
                    elif fname == 'abar':
                        xn_c = np.array(f["abundance"]["xn_c"])
                        a_nuc_rep_c = np.array(f["abundance"]["a_nuc_rep_c"])
                        a_nuc = np.array(f["abundance"]["a_nuc"])
                        a_nuc_tile = np.tile(a_nuc,(xn_c.shape[0],xn_c.shape[1],xn_c.shape[2],1))
                        yn_c = np.empty(xn_c.shape)
                        yn_c[:,:,:,:-1]=xn_c[:,:,:,:-1]/a_nuc_tile[:,:,:,:]
                        yn_c[:,:,:,-1]=xn_c[:,:,:,-1]/a_nuc_rep_c[:,:,:]
                        ytot = np.sum(yn_c,axis=3)
                        atot = np.sum(xn_c,axis=3)
                        abar = np.divide(atot,ytot)
                        dat_1 = abar[:,:,:].transpose()
                    elif fname in ('e_rms_1','e_rms_2','e_rms_3','e_rms_4'):
                        dims = f["mesh"]["array_dimensions"]
                        n_groups = f["radiation"]["raddim"][0]
                        n_species= f["radiation"]["raddim"][1]
                        n_hyperslabs = f["mesh"]["nz_hyperslabs"][()]
                        energy_edge = f["radiation"]["unubi"][()]
                        energy_center = f["radiation"]["unui"][()]
                        d_energy=[]
                        for i in range(0,n_groups):
                            d_energy.append(energy_edge[i+1]-energy_edge[i])
                        d_energy=np.array(d_energy)
                        e3de = energy_center**3*d_energy
                        e5de = energy_center**5*d_energy

                        psi0_c = f["radiation"]["psi0_c"][:]
                        row = np.empty((n_species,int(dims[2]/n_hyperslabs),dims[1],dims[0]))
                        for n in range(0,n_species):
                            numerator=np.sum(psi0_c[:,:,:,n]*e5de,axis=3)
                            denominator=np.sum(psi0_c[:,:,:,n]*e3de,axis=3)
                            row[n][:][:][:]=np.sqrt(numerator/(denominator+1e-100))
                            species = int(fname[-1])-1
                            dat_1 = row[species,:,:,:].transpose()
                    elif fname in ('lumin_1','lumin_2','lumin_3','lumin_4',
                                   'num_lumin_1','num_lumin_2','num_lumin_3','num_lumin_4'):
                        dims = f["mesh"]["array_dimensions"]
                        n_groups = f["radiation"]["raddim"][0]
                        n_hyperslabs = f["mesh"]["nz_hyperslabs"][()]
                        ergmev = 1.602177e-6
                        cvel=2.99792458e10
                        h = 4.13567e-21
                        ecoef = 4.0 * np.pi * ergmev/(h*cvel)**3
                        radius = f['mesh']['x_ef'][()]
                        agr_e = f['fluid']['agr_e'][()]
                        cell_area_GRcorrected = 4*np.pi*radius**2/agr_e**4
                        psi1_e = f['radiation']['psi1_e']
                        energy_edge = f["radiation"]["unubi"][()]
                        energy_center = f["radiation"]["unui"][()]
                        d_energy=[]
                        for i in range(0,n_groups):
                            d_energy.append(energy_edge[i+1]-energy_edge[i])
                        d_energy=np.array(d_energy)
                        e2de = energy_center**2*d_energy
                        e3de = energy_center**3*d_energy
                        e5de = energy_center**5*d_energy
                        species = int(fname[-1])-1
                        if fname in ('lumin_1','lumin_2','lumin_3','lumin_4'):
                            lumin = np.sum(psi1_e[:,:,:,species]*e3de, axis=3)*np.tile(cell_area_GRcorrected[1:dims[0]+1],(int(dims[2]/n_hyperslabs),dims[1],1))*(cvel*ecoef*1e-51)
                        elif fname in ('num_lumin_1','num_lumin_2','num_lumin_3','num_lumin_4'):
                            lumin = np.sum(psi1_e[:,:,:,species]*e2de, axis=3)*np.tile(cell_area_GRcorrected[1:dims[0]+1],(int(dims[2]/n_hyperslabs),dims[1],1))*(cvel*ecoef*1e-51)
                        dat_1 = lumin[:,:,:].transpose()
                        
                        
                    elif fname == 'shock':
                        sys.exit('Error: Currently disfuctional, please select another option')
                        p = np.array(f["fluid"]["press"][:])
                        p = p.reshape(int(f["fluid"]["press"][:][0].size/f["fluid"]["press"][0][0].size), f["fluid"]["press"][0][0].size)
                        u = np.array(f["fluid"]["u_c"][:])
                        u = u.reshape(int(f["fluid"]["u_c"][:][0].size/f["fluid"]["u_c"][0][0].size), f["fluid"]["u_c"][0][0].size)
                        dim1 = f["fluid"]["u_c"][:][0].size/f["fluid"]["u_c"][0][0].size
                        dim2 = f["fluid"]["u_c"][0][0].size
                        f = np.zeros([dim1+1,dim2+1])
                        small = 1.0e-54
                        
                        
                        def shock(p, u, n):
                            pf = np.zeros(dim2)
                            final = np.zeros(dim2)
                            delp1 = np.zeros(dim2)
                            delp2 = np.zeros(dim2)
                            for i in range(dim2-2):
                                if i >= 2 and i < dim2-2:
                                    delp1[i] = float(p[n,i+1] - p[n,i-1])
                                    delp2[i] = float(p[n,i+2] - p[n,i-2])
                                    shk = delp1[i]/(delp2[i]+small)
                                    pf[i] = max(0.0, min(1.0, 10*(shk-0.75)))
                                    exp1 = abs(delp1[i])/(min(p[n,i+1],p[n,i-1])+small)
                                    if exp1 < (1/3):
                                        pf[i] = 0.0
                                    if u[n,i+1]-u[n,i-1]>0:
                                        pf[i] = 0.0
                            for i in range(dim2-2):
                                if i < dim2-2 and delp1[i]<0:
                                    final[i] = max(pf[i], pf[i-1])
                                elif i < dim2-2:
                                    final[i] = max(pf[i],pf[i-1])
                                if i == dim2-3 or i == dim2-4:
                                    if final[i-1]==0 and final[i]:
                                        final[i]=0
                            f[n,:dim2] = final
                            f = f.reshape([dims[2],dims[1]+1,dims[0]+1])
                        
                        for n in range(dim1):
                            shock(p,u,n)
                                
                                
                                
                                                                                                
                                                                                                        
                                                                                                        
                                                                                                        
                    if (f["mesh"]["array_dimensions"][2] > 1):
                        data = dat_1[:-2,:,:] #Clips off ghost zones for 3D
                    else:
                        data = dat_1[:,:,:]
                        data = np.ma.masked_where(mask == 0.0 , data) #Masks
                        data = np.ma.filled(data, fill_value = 0.0) #Replaces masked value with 0
                        
                        
                        
                    yield field, chunk.objs[0], data
        pass

    def _read_chunk_data(self, chunk, fields):
        # This reads the data from a single chunk without doing any selection,
        # and is only used for caching data that might be used by multiple
        # different selectors later. For instance, this can speed up ghost zone
        # computation.
        pass

