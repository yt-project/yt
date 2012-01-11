"""
Skeleton objects that represent a few fundamental yt data types.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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

import abc

class ContainerClass(object):
    pass

class MinimalRepresentation(object):
    __metaclass__ = abc.ABCMeta

    def _update_attrs(self, obj, attr_list):
        for attr in attr_list:
            setattr(self, attr, getattr(obj, attr, None))
        if hasattr(obj, "pf"):
            self.output_hash = obj.pf._hash()

    def __init__(self, obj):
        self._update_attrs(obj, self._attr_list)

    @abc.abstractmethod
    def _generate_post(self):
        pass

    @abc.abstractproperty
    def _attr_list(self):
        pass

    def _return_filtered_object(self, attrs):
        new_attrs = tuple(attr for attr in self._attr_list
                          if attr not in attrs)
        new_class = type('Filtered%s' % self.__class__.__name__,
                         (FilteredRepresentation,),
                         {'_attr_list': new_attrs})
        return new_class(self)

    @property
    def _attrs(self):
        return dict( ((attr, getattr(self, attr)) for attr in self._attr_list) )

    @classmethod
    def _from_metadata(cls, metadata):
        cc = ContainerClass()
        for a, v in metadata.values():
            setattr(cc, a, v)
        return cls(cc)

class FilteredRepresentation(MinimalRepresentation):
    def _generate_post(self):
        raise RuntimeError

class MinimalStaticOutput(MinimalRepresentation):
    _attr_list = ("dimensionality", "refine_by", "domain_dimensions",
                  "current_time", "domain_left_edge", "domain_right_edge",
                  "unique_identifier", "current_redshift", "output_hash",
                  "cosmological_simulation", "omega_matter", "omega_lambda",
                  "hubble_constant", "name")

    def __init__(self, obj):
        super(MinimalStaticOutput, self).__init__(obj)
        self.output_hash = obj._hash()
        self.name = str(obj)

    def _generate_post(self):
        metadata = self._attrs
        chunks = []
        return metadata, chunks

class MinimalMappableData(MinimalRepresentation):

    weight = "None"
    _attr_list = ("field_data", "field", "weight", "axis", "output_hash")

    def _generate_post(self):
        nobj = self._return_filtered_object(("field_data",))
        metadata = nobj._attrs
        chunks = [(arr, self.field_data[arr]) for arr in self.field_data]
        return (metadata, chunks)

class MinimalProjectionData(MinimalMappableData):

    def __init__(self, obj):
        super(MinimalProjectionData, self).__init__(obj)
        self.type = "proj"
