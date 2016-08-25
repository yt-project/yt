"""
OWLS fields




"""
from __future__ import absolute_import

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import numpy as np
from . import owls_ion_tables as oit

from yt.funcs import \
    mylog, download_file
from yt.config import ytcfg
from yt.fields.particle_fields import \
    add_volume_weighted_smoothed_field
from yt.fields.species_fields import \
    add_species_field_by_fraction, \
    add_species_field_by_density
from yt.frontends.sph.fields import \
    SPHFieldInfo



class OWLSFieldInfo(SPHFieldInfo):

    _ions = ("c1", "c2", "c3", "c4", "c5", "c6",
             "fe2", "fe17", "h1", "he1", "he2", "mg1", "mg2", "n2", 
             "n3", "n4", "n5", "n6", "n7", "ne8", "ne9", "ne10", "o1", 
             "o6", "o7", "o8", "si2", "si3", "si4", "si13")

    _elements = ("H", "He", "C", "N", "O", "Ne", "Mg", "Si", "Fe")

    _num_neighbors = 48

    _add_elements = ("PartType0", "PartType4")

    _add_ions = ("PartType0")


    def __init__(self, *args, **kwargs):
        
        new_particle_fields = (
            ("Hydrogen", ("", ["H_fraction"], None)),
            ("Helium", ("", ["He_fraction"], None)),
            ("Carbon", ("", ["C_fraction"], None)),
            ("Nitrogen", ("", ["N_fraction"], None)),
            ("Oxygen", ("", ["O_fraction"], None)),
            ("Neon", ("", ["Ne_fraction"], None)),
            ("Magnesium", ("", ["Mg_fraction"], None)),
            ("Silicon", ("", ["Si_fraction"], None)),
            ("Iron", ("", ["Fe_fraction"], None))
            )

        self.known_particle_fields += new_particle_fields

        super(OWLSFieldInfo,self).__init__( *args, **kwargs )

        # This enables the machinery in yt.fields.species_fields
        self.species_names += list(self._elements)

    def setup_particle_fields(self, ptype):
        """ additional particle fields derived from those in snapshot.
        we also need to add the smoothed fields here b/c setup_fluid_fields
        is called before setup_particle_fields. """ 

        smoothed_suffixes = ("_number_density", "_density", "_mass")

        # we add particle element fields for stars and gas
        #-----------------------------------------------------
        if ptype in self._add_elements:


            # this adds the particle element fields
            # X_density, X_mass, and X_number_density
            # where X is an item of self._elements.
            # X_fraction are defined in snapshot
            #-----------------------------------------------
            for s in self._elements:
                add_species_field_by_fraction(self, ptype, s,
                                              particle_type=True)

        # this needs to be called after the call to 
        # add_species_field_by_fraction for some reason ...
        # not sure why yet. 
        #-------------------------------------------------------
        if ptype == 'PartType0':
            ftype='gas'
        elif ptype == 'PartType1':
            ftype='dm'
        elif ptype == 'PartType2':
            ftype='PartType2'
        elif ptype == 'PartType3':
            ftype='PartType3'
        elif ptype == 'PartType4':
            ftype='star'
        elif ptype == 'PartType5':
            ftype='BH'
        elif ptype == 'all':
            ftype='all'
        
        super(OWLSFieldInfo,self).setup_particle_fields(
            ptype, num_neighbors=self._num_neighbors, ftype=ftype)


        # and now we add the smoothed versions for PartType0
        #-----------------------------------------------------
        if ptype == 'PartType0':

            # we only add ion fields for gas.  this takes some 
            # time as the ion abundances have to be interpolated
            # from cloudy tables (optically thin)
            #-----------------------------------------------------
    

            # this defines the ion density on particles
            # X_density for all items in self._ions
            #-----------------------------------------------
            self.setup_gas_ion_density_particle_fields( ptype )

            # this adds the rest of the ion particle fields
            # X_fraction, X_mass, X_number_density
            #-----------------------------------------------
            for ion in self._ions:

                # construct yt name for ion
                #---------------------------------------------------
                if ion[0:2].isalpha():
                    symbol = ion[0:2].capitalize()
                    roman = int(ion[2:])
                else:
                    symbol = ion[0:1].capitalize()
                    roman = int(ion[1:])

                if (ptype, symbol + "_fraction") not in self.field_aliases:
                    continue

                pstr = "_p" + str(roman-1)
                yt_ion = symbol + pstr

                # add particle field
                #---------------------------------------------------
                add_species_field_by_density(self, ptype, yt_ion,
                                             particle_type=True)


            # add smoothed ion fields
            #-----------------------------------------------
            for ion in self._ions:

                # construct yt name for ion
                #---------------------------------------------------
                if ion[0:2].isalpha():
                    symbol = ion[0:2].capitalize()
                    roman = int(ion[2:])
                else:
                    symbol = ion[0:1].capitalize()
                    roman = int(ion[1:])

                if (ptype, symbol + "_fraction") not in self.field_aliases:
                    continue

                pstr = "_p" + str(roman-1)
                yt_ion = symbol + pstr

                loaded = []
                for sfx in smoothed_suffixes:
                    fname = yt_ion + sfx
                    fn = add_volume_weighted_smoothed_field( 
                        ptype, "particle_position", "particle_mass",
                        "smoothing_length", "density", fname, self,
                        self._num_neighbors)
                    loaded += fn

                    self.alias(("gas", fname), fn[0])

                self._show_field_errors += loaded
                self.find_dependencies(loaded)



    def setup_gas_ion_density_particle_fields( self, ptype ):
        """ Sets up particle fields for gas ion densities. """ 

        # loop over all ions and make fields
        #----------------------------------------------
        for ion in self._ions:

            # construct yt name for ion
            #---------------------------------------------------
            if ion[0:2].isalpha():
                symbol = ion[0:2].capitalize()
                roman = int(ion[2:])
            else:
                symbol = ion[0:1].capitalize()
                roman = int(ion[1:])

            if (ptype, symbol + "_fraction") not in self.field_aliases:
                continue

            pstr = "_p" + str(roman-1)
            yt_ion = symbol + pstr
            ftype = ptype

            # add ion density field for particles
            #---------------------------------------------------
            fname = yt_ion + '_density'
            dens_func = self._create_ion_density_func( ftype, ion )
            self.add_field( (ftype, fname),
                            function = dens_func, 
                            units=self.ds.unit_system["density"],
                            particle_type=True )            
            self._show_field_errors.append( (ftype,fname) )



        
    def _create_ion_density_func( self, ftype, ion ):
        """ returns a function that calculates the ion density of a particle. 
        """ 

        def get_owls_ion_density_field(ion, ftype, itab):
            def _func(field, data):

                # get element symbol from ion string. ion string will 
                # be a member of the tuple _ions (i.e. si13)
                #--------------------------------------------------------
                if ion[0:2].isalpha():
                    symbol = ion[0:2].capitalize()
                else:
                    symbol = ion[0:1].capitalize()

                # mass fraction for the element
                #--------------------------------------------------------
                m_frac = data[ftype, symbol+"_fraction"]

                # get nH and T for lookup
                #--------------------------------------------------------
                log_nH = np.log10( data["PartType0", "H_number_density"] )
                log_T = np.log10( data["PartType0", "Temperature"] )

                # get name of owls_ion_file for given ion
                #--------------------------------------------------------
                itab.set_iz( data.ds.current_redshift )

                # find ion balance using log nH and log T
                #--------------------------------------------------------
                i_frac = itab.interp( log_nH, log_T )
                return data[ftype,"Density"] * m_frac * i_frac 
            return _func
            
        ion_path = self._get_owls_ion_data_dir()
        fname = os.path.join( ion_path, ion+".hdf5" )
        itab = oit.IonTableOWLS( fname )
        return get_owls_ion_density_field(ion, ftype, itab)





    # this function sets up the X_mass, X_density, X_fraction, and
    # X_number_density fields where X is the name of an OWLS element.
    #-------------------------------------------------------------
    def setup_fluid_fields(self):

        return



    # this function returns the owls_ion_data directory. if it doesn't
    # exist it will download the data from http://yt-project.org/data
    #-------------------------------------------------------------
    def _get_owls_ion_data_dir(self):

        txt = "Attempting to download ~ 30 Mb of owls ion data from %s to %s."
        data_file = "owls_ion_data.tar.gz"
        data_url = "http://yt-project.org/data"

        # get test_data_dir from yt config (ytcgf)
        #----------------------------------------------
        tdir = ytcfg.get("yt","test_data_dir")

        # set download destination to tdir or ./ if tdir isnt defined
        #----------------------------------------------
        if tdir == "/does/not/exist":
            data_dir = "./"
        else:
            data_dir = tdir            


        # check for owls_ion_data directory in data_dir
        # if not there download the tarball and untar it
        #----------------------------------------------
        owls_ion_path = os.path.join( data_dir, "owls_ion_data" )

        if not os.path.exists(owls_ion_path):
            mylog.info(txt % (data_url, data_dir))
            fname = data_dir + "/" + data_file
            download_file(os.path.join(data_url, data_file), fname)

            cmnd = "cd " + data_dir + "; " + "tar xf " + data_file
            os.system(cmnd)


        if not os.path.exists(owls_ion_path):
            raise RuntimeError("Failed to download owls ion data.")

        return owls_ion_path
