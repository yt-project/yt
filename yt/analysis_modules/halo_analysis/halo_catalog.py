"""
HaloCatalog object



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import h5py
import numpy as np
import os

from yt.funcs import \
     ensure_dir, \
     mylog
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, \
    parallel_blocking_call, \
    parallel_root_only, \
    parallel_objects
     
from .halo_object import \
     Halo
from .operator_registry import \
    callback_registry, \
    quantity_registry, \
    hf_registry

class HaloCatalog(ParallelAnalysisInterface):
    r"""Create a HaloCatalog: an object that allows for the creation and association
    of data with a set of halo objects.

    Parameters
    ----------
    halos_pf : str
        Dataset created by a halo finder.  If None, a halo finder should be 
        provided with the finder_method keyword.
    data_pf : str
        Dataset created by a simulation.
    data_source : data container
        Data container associated with either the halos_pf or the data_pf.
    finder_method : str
        Halo finder to be used if no halos_pf is given.
    output_dir : str
        The top level directory into which analysis output will be written.
        Default: "."

    Examples
    --------
    
    """
    
    def __init__(self, halos_pf=None, data_pf=None, 
                 data_source=None, finder_method=None, 
                 output_dir="."):
        ParallelAnalysisInterface.__init__(self)
        self.halos_pf = halos_pf
        self.data_pf = data_pf
        self.finder_method = finder_method
        self.output_dir = ensure_dir(output_dir)

        if halos_pf is None:
            if data_pf is None:
                raise RuntimeError("Must specify a halos_pf, data_pf, or both.")
            if finder_method is None:
                raise RuntimeError("Must specify a halos_pf or a finder_method.")
        if data_source is None:
            if halos_pf is not None:
                data_source = halos_pf.h.all_data()
            else:
                data_source = data_pf.h.all_data()
        self.data_source = data_source

        self.add_default_quantities()
        self.callbacks = []

    def add_callback(self, callback, *args, **kwargs):
        callback = callback_registry.find(callback, *args, **kwargs)
        if "output_dir" in kwargs is not None:
            ensure_dir(os.path.join(self.output_dir, kwargs["output_dir"]))
        self.callbacks.append(callback)

    def add_quantity(self, key, field_type=None, *args, **kwargs):
        if field_type is None:
            quantity = quantity_registry.find(key, *args, **kwargs)
        elif (field_type, key) in self.halos_pf.field_info:
            quantity = (field_type, key)
        else:
            raise RuntimeError("HaloCatalog quantity must be a registered function or a field of a known type.")
        self.quantities.append((key, quantity))

    def add_filter(self, filter, *args, **kwargs):
        filter = callback_registry.find(filter, *args, **kwargs)
        self.callbacks.append(filter)

    def run(self, njobs=-1, filename=None):
        self.halo_list = []

        if self.halos_pf is None:
            # this is where we would do halo finding and assign halos_pf to 
            # the dataset that we have just created.
            raise NotImplementedError

        self.n_halos = self.data_source["particle_identifier"].size
        for i in parallel_objects(xrange(self.n_halos), njobs=njobs):
            new_halo = Halo(self)
            for key, quantity in self.quantities:
                if quantity in self.halos_pf.field_info:
                    new_halo.quantities[key] = self.data_source[quantity][i]
                elif callable(quantity):
                    new_halo.quantities[key] = quantity(new_halo)

            for callback in self.callbacks:
                callback(new_halo)

            self.halo_list.append(new_halo.quantities)
            
            del new_halo

        self.halo_list = self.comm.par_combine_object(self.halo_list,
                                                      datatype="list",
                                                      op="cat")
        self.halo_list.sort(key=lambda a:a['particle_identifier'].to_ndarray())
            
        if filename is not None:
            self.save_catalog(filename)

    @parallel_root_only
    def save_catalog(self, filename):
        "Write out hdf5 file with all halo quantities."

        mylog.info("Saving halo catalog to %s." % filename)
        out_file = h5py.File(filename, 'w')
        field_data = np.empty(self.n_halos)
        for key, quantity in self.quantities:
            for i in xrange(self.n_halos):
                field_data[i] = self.halo_list[i][key]
            out_file.create_dataset(key, data=field_data)
        out_file.close()

    def add_default_quantities(self):
        self.quantities = []
        self.add_quantity("particle_identifier", field_type="halos")
        self.add_quantity("particle_mass", field_type="halos")
        self.add_quantity("particle_position_x", field_type="halos")
        self.add_quantity("particle_position_y", field_type="halos")
        self.add_quantity("particle_position_z", field_type="halos")
        self.add_quantity("virial_radius", field_type="halos")
        
