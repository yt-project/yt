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
     parallel_objects
     
from .halo_object import \
     Halo
from .operator_registry import \
     callback_registry, \
     filter_registry, \
     hf_registry, \
     quantity_registry

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
                 output_dir="halo_catalogs/catalog"):
        ParallelAnalysisInterface.__init__(self)
        self.halos_pf = halos_pf
        self.data_pf = data_pf
        self.finder_method = finder_method
        self.output_dir = ensure_dir(output_dir)
        if os.path.basename(self.output_dir) != ".":
            self.output_prefix = os.path.basename(self.output_dir)
        else:
            self.output_prefix = "catalog"

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

        # all of the analysis actions to be performed: callbacks, filters, and quantities
        self.actions = []
        # fields to be written to the halo catalog
        self.quantities = []
        self.add_default_quantities()

    def add_callback(self, callback, *args, **kwargs):
        callback = callback_registry.find(callback, *args, **kwargs)
        if "output_dir" in kwargs is not None:
            ensure_dir(os.path.join(self.output_dir, kwargs["output_dir"]))
        self.actions.append(("callback", callback))

    def add_quantity(self, key, *args, **kwargs):
        if "field_type" in kwargs:
            field_type = kwargs.pop("field_type")
        else:
            field_type = None
        if field_type is None:
            quantity = quantity_registry.find(key, *args, **kwargs)
        elif (field_type, key) in self.halos_pf.field_info:
            quantity = (field_type, key)
        else:
            raise RuntimeError("HaloCatalog quantity must be a registered function or a field of a known type.")
        self.quantities.append(key)
        self.actions.append(("quantity", (key, quantity)))

    def add_filter(self, halo_filter, *args, **kwargs):
        halo_filter = filter_registry.find(halo_filter, *args, **kwargs)
        self.actions.append(("filter", halo_filter))

    def create(self, save_halos=False, save_catalog=True, njobs=-1, dynamic=False):
        self._run(save_halos, save_catalog, njobs=njobs, dynamic=dynamic)

    def load(self, save_halos=True, save_catalog=False, njobs=-1, dynamic=False):
        self._run(save_halos, save_catalog, njobs=njobs, dynamic=dynamic)
        
    @parallel_blocking_call
    def _run(self, save_halos, save_catalog, njobs=-1, dynamic=False):
        self.catalog = []
        if save_halos: self.halo_list = []

        if self.halos_pf is None:
            # this is where we would do halo finding and assign halos_pf to 
            # the dataset that we have just created.
            raise NotImplementedError

        my_index = np.argsort(self.data_source["particle_identifier"])
        for i in parallel_objects(my_index, njobs=njobs, dynamic=dynamic):
            new_halo = Halo(self)
            halo_filter = True
            for action_type, action in self.actions:
                if action_type == "callback":
                    action(new_halo)
                elif action_type == "filter":
                    halo_filter = action(new_halo)
                    if not halo_filter: break
                elif action_type == "quantity":
                    key, quantity = action
                    if quantity in self.halos_pf.field_info:
                        new_halo.quantities[key] = \
                          self.data_source[quantity][int(i)].in_cgs()
                    elif callable(quantity):
                        new_halo.quantities[key] = quantity(new_halo)
                else:
                    raise RuntimeError("Action must be a callback, filter, or quantity.")

            if halo_filter:
                self.catalog.append(new_halo.quantities)

            if save_halos and halo_filter:
                self.halo_list.append(new_halo)
            else:
                del new_halo

        self.catalog.sort(key=lambda a:a['particle_identifier'].to_ndarray())
        if save_catalog:
            self.save_catalog()

    def save_catalog(self):
        "Write out hdf5 file with all halo quantities."

        filename = os.path.join(self.output_dir, "%s.%d.h5" %
                                (self.output_prefix, self.comm.rank))
        n_halos = len(self.catalog)
        mylog.info("Saving halo catalog (%d halos) to %s." %
                   (n_halos, os.path.join(self.output_dir, 
                                         self.output_prefix)))
        out_file = h5py.File(filename, 'w')
        for attr in ["current_redshift", "current_time",
                     "domain_dimensions",
                     "cosmological_simulation", "omega_lambda",
                     "omega_matter", "hubble_constant"]:
            out_file.attrs[attr] = getattr(self.halos_pf, attr)
        for attr in ["domain_left_edge", "domain_right_edge"]:
            out_file.attrs[attr] = getattr(self.halos_pf, attr).in_cgs()
        out_file.attrs["data_type"] = "halo_catalog"
        out_file.attrs["num_halos"] = n_halos
        if n_halos > 0:
            field_data = np.empty(n_halos)
            for key in self.quantities:
                units = ""
                if hasattr(self.catalog[0][key], "units"):
                    units = str(self.catalog[0][key].units)
                for i in xrange(n_halos):
                    field_data[i] = self.catalog[i][key]
                dataset = out_file.create_dataset(str(key), data=field_data)
                dataset.attrs["units"] = units
        out_file.close()

    def add_default_quantities(self):
        self.add_quantity("particle_identifier", field_type="halos")
        self.add_quantity("particle_mass", field_type="halos")
        self.add_quantity("particle_position_x", field_type="halos")
        self.add_quantity("particle_position_y", field_type="halos")
        self.add_quantity("particle_position_z", field_type="halos")
        self.add_quantity("virial_radius", field_type="halos")
        
