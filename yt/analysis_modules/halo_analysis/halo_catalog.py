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

class HaloCatalog(object):
    def __init__(self, pf, finding_method, data_source = None):
        self.pf = pf
        if not callable(finding_method):
            finding_method = hf_registry.find(finding_method)
        self.finding_method = finding_method

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
