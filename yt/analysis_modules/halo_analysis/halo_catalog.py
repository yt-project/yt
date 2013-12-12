"""
Halo Catalog object



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from .operator_registry import \
    callback_registry, \
    quantity_registry, \
    hf_registry

class HaloCatalog(object):
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

    Examples
    --------
    
    """
    
    def __init__(self, halos_pf=None, data_pf=None, data_source=None, finder_method=None):
        self.halos_pf = halos_pf
        self.data_pf = data_pf
        self.finder_method = finder_method

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

        self.values = []
        self.callbacks = []

    def add_callback(self, callback, *args, **kwargs):
        callback = callback_registry.find(callback, *args, **kwargs)
        self.callbacks.append(callback)

    def add_quantity(self, quantity, *args, **kwargs):
        quantity = quantity_registry.find(quantity, *args, **kwargs)
        self.callbacks.append(quantity)

    def add_filter(self, filter, *args, **kwargs):
        filter = callback_registry.find(filter, *args, **kwargs)
        self.callbacks.append(filter)

    def run(self):

        if halos_pf is None:
            # this is where we would do halo finding and assign halos_pf to 
            # the dataset that we have just created.
            raise NotImplementedError
        
        self.run_callbacks(halo_list)

    def run_callbacks(self, halo_list):
        for cb in self.callbacks:
            cb.initialize(self)
        for halo in halo_list:
            if all(cb(self, halo) for cb in self.callbacks):
                self.values.append(halo.quantities)
        for cb in self.callbacks:
            cb.finalize(self)

