"""
Various definitions for various other modules and routines - modified
from Orion frontend.


Authors:
 * J. S. Oishi 
 * Chris Malone 


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from yt.funcs import *

parameterTypes = {"dimensionality": int,
                 "refine_by": int
                 }

# converts the Maestro inputs file name to the Enzo/yt name expected
# throughout the code. key is Maestro name, value is Enzo/yt equivalent
maestro2enzoDict = {"dm_in": "dimensionality",
                    "ref_ratio": "refine_by"
                  }

yt2maestroFieldsDict = {}
maestro2ytFieldsDict = {}

maestro_FAB_header_pattern = r"^FAB \(\((\d+), \([0-9 ]+\)\),\(\d+, \(([0-9 ]+)\)\)\)\(\((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\)\) (\d+)\n"
