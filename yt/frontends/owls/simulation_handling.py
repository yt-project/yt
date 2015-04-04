"""
OWLSSimulation class and member functions.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013-2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os

from yt.frontends.gadget.simulation_handling import \
    GadgetSimulation

class OWLSSimulation(GadgetSimulation):
    r"""
    Initialize an OWLS Simulation object.

    Upon creation, the parameter file is parsed and the time and redshift
    are calculated and stored in all_outputs.  A time units dictionary is
    instantiated to allow for time outputs to be requested with physical
    time units.  The get_time_series can be used to generate a
    DatasetSeries object.

    parameter_filename : str
        The simulation parameter file.
    find_outputs : bool
        If True, the OutputDir directory is searched for datasets.  
        Time and redshift information are gathered by temporarily 
        instantiating each dataset.  This can be used when simulation 
        data was created in a non-standard way, making it difficult 
        to guess the corresponding time and redshift information.
        Default: False.

    Examples
    --------
    >>> import yt
    >>> es = yt.simulation("my_simulation.par", "OWLS")
    >>> es.get_time_series()
    >>> for ds in es:
    ...     print ds.current_time

    """

    def __init__(self, parameter_filename, find_outputs=False):
        GadgetSimulation.__init__(self, parameter_filename,
                                  find_outputs=find_outputs)

    def _snapshot_format(self, index=None):
        """
        The snapshot filename for a given index.  Modify this for different 
        naming conventions.
        """

        if self.parameters["OutputDir"].startswith("/"):
            data_dir = self.parameters["OutputDir"]
        else:
            data_dir = os.path.join(self.directory,
                                    self.parameters["OutputDir"])
        if self.parameters["NumFilesPerSnapshot"] > 1:
            suffix = ".0"
        else:
            suffix = ""
        if self.parameters["SnapFormat"] == 3:
            suffix += ".hdf5"
        if index is None:
            count = "*"
        else:
            count = "%03d" % index
        keyword = "%s_%s" % (self.parameters["SnapshotFileBase"], count)
        filename = os.path.join(keyword, "%s%s" % (keyword, suffix))
        return os.path.join(data_dir, filename)
