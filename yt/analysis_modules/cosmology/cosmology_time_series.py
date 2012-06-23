"""
CosmologyTimeSeries class and member functions.

Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: Michigan State University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2012 Britton Smith.  All Rights Reserved.

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

from yt.convenience import \
    simulation
from yt.funcs import *
from yt.data_objects.time_series import \
    TimeSeriesData
from yt.utilities.parameter_file_storage import \
    simulation_time_series_registry

class CosmologyTimeSeries(TimeSeriesData):
    """
    Class for splicing together datasets to extend over a
    cosmological distance.
    """

    def __init__(self, parameter_filename, simulation_type):
        if simulation_type not in simulation_time_series_registry:
            raise YTSimulationNotIdentified(simulation_type)

        self.parameter_filename = parameter_filename
        self.simulation_type = simulation_type
