"""
OWLS-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import numpy as np
import owls_ion_tables as oit

from yt.funcs import *

from yt.fields.field_info_container import \
    FieldInfoContainer

from .definitions import \
    gadget_ptypes, \
    ghdf5_ptypes

from yt.config import ytcfg
from yt.utilities.physical_constants import mh
from yt.fields.species_fields import \
    add_species_field_by_fraction, \
    add_species_field_by_density

from yt.fields.particle_fields import \
    add_volume_weighted_smoothed_field


# Here are helper functions for things like vector fields and so on.

def _get_conv(cf):
    def _convert(data):
        return data.convert(cf)
    return _convert

class SPHFieldInfo(FieldInfoContainer):
    known_other_fields = ()

    known_particle_fields = (
        ("Mass", ("code_mass", ["particle_mass"], None)),
        ("Masses", ("code_mass", ["particle_mass"], None)),
        ("Coordinates", ("code_length", ["particle_position"], None)),
        ("Velocity", ("code_velocity", ["particle_velocity"], None)),
        ("Velocities", ("code_velocity", ["particle_velocity"], None)),
        ("ParticleIDs", ("", ["particle_index"], None)),
        ("InternalEnergy", ("", ["thermal_energy"], None)),
        ("SmoothingLength", ("code_length", ["smoothing_length"], None)),
        ("Density", ("code_mass / code_length**3", ["density"], None)),
        ("MaximumTemperature", ("K", [], None)),
        ("Temperature", ("K", ["temperature"], None)),
        ("Epsilon", ("code_length", [], None)),
        ("Metals", ("code_metallicity", ["metallicity"], None)),
        ("Phi", ("code_length", [], None)),
        ("FormationTime", ("code_time", ["creation_time"], None)),
    )




class OWLSFieldInfo(SPHFieldInfo):

    _ions = ("c1", "c2", "c3", "c4", "c5", "c6",
             "fe2", "fe17", "h1", "he1", "he2", "mg1", "mg2", "n2", 
             "n3", "n4", "n5", "n6", "n7", "ne8", "ne9", "ne10", "o1", 
             "o6", "o7", "o8", "si2", "si3", "si4", "si13")

    _elements = ("H", "He", "C", "N", "O", "Ne", "Mg", "Si", "Fe")

    # override
    #--------------------------------------------------------------
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






    def setup_particle_fields(self, ptype):
        """ additional particle fields that are not on disk """ 

        super(OWLSFieldInfo,self).setup_particle_fields(ptype)

        if ptype == "PartType0":

            # this adds the particle element fields
            # X_density, X_mass, and X_number_density
            # where X is an element of self._elements.
            # X_fraction are defined in file
            #-----------------------------------------------
            for s in self._elements:
                add_species_field_by_fraction(self, ptype, s,
                                              particle_type=True)

            # this defines the ion density on particles
            # X_density 
            #-----------------------------------------------
            self.setup_ion_density_particle_fields()


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

                pstr = "_p" + str(roman-1)
                yt_ion = symbol + pstr

                add_species_field_by_density(self, ptype, yt_ion,
                                             particle_type=True)



            ptype = 'PartType0'
            num_neighbors = 48
            loaded = []

            # add smoothed element fields
            #-----------------------------------------------
            for s in self._elements:
                fname = s + '_number_density'
                fn = add_volume_weighted_smoothed_field( 
                    ptype, "particle_position", "particle_mass",
                    "smoothing_length", "density", fname, self,
                    num_neighbors)
                loaded += fn

                self.alias(("gas", fname), fn[0])

            self._show_field_errors += loaded


            # add smoothed ion fields
            #-----------------------------------------------
            for ion in self._ions:

                if ion[0:2].isalpha():
                    symbol = ion[0:2].capitalize()
                    roman = int(ion[2:])
                else:
                    symbol = ion[0:1].capitalize()
                    roman = int(ion[1:])

                pstr = "_p" + str(roman-1)
                yt_ion = symbol + pstr

                fname = yt_ion + '_number_density'
                fn = add_volume_weighted_smoothed_field( 
                    ptype, "particle_position", "particle_mass",
                    "smoothing_length", "density", fname, self,
                    num_neighbors)
                loaded += fn

                self.alias(("gas", fname), fn[0])

            self._show_field_errors += loaded


            # find dependencies
            #-----------------------------------------------
            self.find_dependencies(loaded)


            return






    def setup_ion_density_particle_fields( self ):
        """ Sets up particle fields for ion densities. """ 

        # loop over all ions and make fields
        #----------------------------------------------
#        for ion in ("h1",): # self._ions:
        for ion in self._ions:

            # construct yt name for ion
            #---------------------------------------------------
            if ion[0:2].isalpha():
                symbol = ion[0:2].capitalize()
                roman = int(ion[2:])
            else:
                symbol = ion[0:1].capitalize()
                roman = int(ion[1:])

            pstr = "_p" + str(roman-1)
            yt_ion = symbol + pstr
            ftype = "PartType0"

            # add ion density field for particles
            #---------------------------------------------------
            fname = yt_ion + '_density'
            dens_func = self._create_ion_density_func( ftype, ion )
            self.add_field( (ftype, fname),
                            function = dens_func, 
                            units="g/cm**3",
                            particle_type=True )            
            self._show_field_errors.append( (ftype,fname) )



        
    def _create_ion_density_func( self, ftype, ion ):
        """ returns a function that calculates the ion density of a particle. 
        """ 

        def _ion_density(field, data):

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
            owls_ion_path = self._get_owls_ion_data_dir()
            fname = os.path.join( owls_ion_path, ion+".hdf5" )

            # create ionization table for this redshift
            #--------------------------------------------------------
            itab = oit.IonTableOWLS( fname )
            itab.set_iz( data.pf.current_redshift )

            # find ion balance using log nH and log T
            #--------------------------------------------------------
            i_frac = itab.interp( log_nH, log_T )
            return data[ftype,"Density"] * m_frac * i_frac 
        
        return _ion_density





    # this function sets up the X_mass, X_density, X_fraction, and
    # X_number_density fields where X is the name of an OWLS element.
    #-------------------------------------------------------------
    def setup_fluid_fields(self):

        return


        # these add the smoothed element fields
        #-----------------------------------------------
        for s in self._elements:
            add_species_field_by_fraction(self, "gas", s)

        

        # this should add the smoothed ion fields
        #-----------------------------------------------
        for ion in self._ions:

            if ion[0:2].isalpha():
                symbol = ion[0:2].capitalize()
                roman = int(ion[2:])
            else:
                symbol = ion[0:1].capitalize()
                roman = int(ion[1:])

            pstr = "_p" + str(roman-1)
            yt_ion = symbol + pstr

            add_species_field_by_density(self, "gas", yt_ion)


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
            fn = download_file(os.path.join(data_url, data_file), fname)

            cmnd = "cd " + data_dir + "; " + "tar xf " + data_file
            os.system(cmnd)


        if not os.path.exists(owls_ion_path):
            raise RuntimeError, "Failed to download owls ion data."

        return owls_ion_path
