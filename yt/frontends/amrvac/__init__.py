"""
API for yt.frontends.amrvac



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _f90nml as f90nml
from yt.extern.six import string_types

def read_amrvac_namelist(parfiles):
    """Read one or more parfiles, and return a unified f90nml.Namelist object.

    This function replicates the patching logic of MPI-AMRVAC where redundant parameters
    only retain last-in-line values EXCEPT `&filelist:base_filename`, which is accumulated.
    When passed a single file, this function acts as a mere wrapper of f90nml.read().

    Parameters
    ----------
    parfiles : str or list
        A file path, or a list of file paths to MPI-AMRVAC configuration parfiles.

    Returns
    -------
    unified_namelist : f90nml.Namelist
        A single namelist object. The class inherits from ordereddict.

    """
    # typechecking
    if isinstance(parfiles, string_types):
        parfiles = [parfiles]
    assert all([isinstance(pf, string_types) for pf in parfiles])

    # first merge the namelists
    namelists = [f90nml.read(parfile) for parfile in parfiles]
    unified_namelist = f90nml.Namelist()
    for nml in namelists:
        unified_namelist.patch(nml)

    # accumulate `&filelist:base_filename`
    base_filename = "".join([nml.get("filelist", {}).get("base_filename", "") for nml in namelists])
    unified_namelist["filelist"]["base_filename"] = base_filename

    return unified_namelist
