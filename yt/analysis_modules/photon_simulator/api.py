"""
API for yt.analysis_modules.photon_simulator.
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from .photon_simulator import \
     PhotonList, \
     EventList

from .photon_models import \
     XSpecThermalModel, \
     XSpecAbsorbModel, \
     TableApecModel, \
     TableAbsorbModel
