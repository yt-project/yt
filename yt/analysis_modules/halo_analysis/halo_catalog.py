"""
Halo Catalog object

Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: MSU
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2013 Britton Smith, Matthew Turk.  All Rights Reserved.

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

from .operator_registry impoort \
    hf_registry, \
    callback_registry, \
    filter_registry, \
    quantity_registry

class HaloCatalog(object):
    def __init__(self, pf, finding_method, data_source = None):
        self.pf = pf
        if not callable(finding_method):
            finding_method = hf_registry.find(finding_method)
        self.finding_method = finding_method
        if data_source is None:
            data_source = pf.h.all_data()
        self.data_source = data_source

        # Now we initialize our callbacks, quantities and filters
        self.filters = []
        self.quantity_operators = []
        self.callbacks = []

    def add_callback(self, callback):
        if not callable(callback):
            callback = callback_registry.find(callback)
        self.callbacks.append(callback)

    def add_filter(self, filter):
        if not callable(filter):
            filter = filter_registry.find(filter)
        self.filters.append(filter)

    def add_quantity(self, quantity):
        if not callable(quantity):
            quantity = quantity_registry.find(quantity)
        self.quantity_operators.append(quantity)

    def run(self):
        # Here's the basic rundown.
        # First we call the halo finding operation.  This is going to be handed
        # the data source, but we assume it already has all its arguments
        # necessary for the finding operation.
        halo_list = self.finding_method(self)
        # This will be our first pass list of halos.
        # Now we filter our halos.  So for each filter and for each halo, we
        # will filter and construct a new set of halos.
        halo_list = self.filter_halos(halo_list)
        # Once we've filtered, we can go ahead and calculate properties.
        self.compute_properties(halo_list)
        # Now all the halos have their properties.
        self.halo_list = halo_list
        self.run_callbacks(halo_list)

    def filter_halos(self, halo_list):
        new_list = []
        for hid in range(len(halo_list)):
            halos = [halo_list.pop(0)]
            for filter in self.filters:
                new_halos = halos
                for halo in new_halos:
                    new_halos += filter(self, halo)
                halos = new_halos
            new_list += halos
        return new_list

    def compute_properties(self, halo_list):
        for halo in halo_list:
            for quantity in self.quantity_operators:
                quantity(self, halo)

    def run_callbacks(self, halo_list):
        for halo in halo_list:
            for callback in self.callbacks:
                callback(self, halo)
