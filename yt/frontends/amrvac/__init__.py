import os

from yt.funcs import ensure_list
from yt.utilities.on_demand_imports import _f90nml as f90nml


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
    parfiles = [os.path.expanduser(pf) for pf in ensure_list(parfiles)]

    # first merge the namelists
    namelists = [f90nml.read(parfile) for parfile in parfiles]
    unified_namelist = f90nml.Namelist()
    for nml in namelists:
        unified_namelist.patch(nml)

    # accumulate `&filelist:base_filename`
    base_filename = "".join(
        [nml.get("filelist", {}).get("base_filename", "") for nml in namelists]
    )
    unified_namelist["filelist"]["base_filename"] = base_filename

    return unified_namelist
