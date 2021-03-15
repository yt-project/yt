import os
from copy import deepcopy
from pathlib import Path

from yt.frontends.amrvac.api import read_amrvac_namelist
from yt.testing import requires_module
from yt.utilities.on_demand_imports import _f90nml as f90nml

test_dir = os.path.dirname(os.path.abspath(__file__))
blast_wave_parfile = os.path.join(test_dir, "sample_parfiles", "bw_3d.par")
modifier_parfile = os.path.join(test_dir, "sample_parfiles", "tvdlf_scheme.par")


@requires_module("f90nml")
def test_read_pathlike():
    read_amrvac_namelist(Path(blast_wave_parfile))


@requires_module("f90nml")
def test_read_one_file():
    """when provided a single file, the function should merely act
    as a wrapper for f90nml.read()"""
    namelist1 = read_amrvac_namelist(blast_wave_parfile)
    namelist2 = f90nml.read(blast_wave_parfile)
    assert namelist1 == namelist2


@requires_module("f90nml")
def test_accumulate_basename():
    """When two (or more) parfiles are passed,
    the filelist:base_filename should be special-cased"""
    namelist_base = f90nml.read(blast_wave_parfile)
    namelist_update = f90nml.read(modifier_parfile)

    namelist_tot1 = read_amrvac_namelist([blast_wave_parfile, modifier_parfile])
    namelist_tot2 = deepcopy(namelist_base)
    namelist_tot2.patch(namelist_update)

    # remove and store the special-case value
    name1 = namelist_tot1["filelist"].pop("base_filename")
    name2 = namelist_tot2["filelist"].pop("base_filename")

    assert (
        name1
        == namelist_base["filelist"]["base_filename"]
        + namelist_update["filelist"]["base_filename"]
    )
    assert name2 == namelist_update["filelist"]["base_filename"]
    assert name1 != name2

    # test equality for the rest of the namelist
    assert namelist_tot1 == namelist_tot2
