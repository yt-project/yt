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

from yt.funcs import issue_deprecation_warning

issue_deprecation_warning("The photon_simulator module is deprecated. Please use pyXSIM "
                          "(http://hea-www.cfa.harvard.edu/~jzuhone/pyxsim) instead.")

from .photon_models import \
     PhotonModel, \
     ThermalPhotonModel

from .photon_simulator import \
     PhotonList, \
     EventList, \
     merge_files, \
     convert_old_file

from .spectral_models import \
     SpectralModel, \
     XSpecThermalModel, \
     XSpecAbsorbModel, \
     TableApecModel, \
     TableAbsorbModel
