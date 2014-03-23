"""
Import the components of the volume rendering extension



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.orientation import Orientation
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import ParallelAnalysisInterface


class Camera(ParallelAnalysisInterface):

    r"""    """


    def __init__(self):
        """Initialize a Camera Instance"""
        ParallelAnalysisInterface.__init__(self)
        self.position = None
        self.focus = None
        self.orientation = None
        self.light = None

