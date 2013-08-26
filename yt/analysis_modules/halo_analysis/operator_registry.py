"""
Operation registry class

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

class OperatorRegistry(dict):
    def find(self, op):
        args = None
        if isinstance(op, tuple):
            op, args = op
        if not callable(op):
            # Lookup, assuming string or hashable object
            op = self[op]
        # We assume that if we are fed args, and a callable, it is arguments
        # for a class instantiation.
        if args is not None:
            op = op(*args)
        # This should at this point be a final operation.
        return op

callback_registry = OperatorRegistry()
hf_registry = OperatorRegistry()
