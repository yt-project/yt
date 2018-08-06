"""
light cone generator test



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import \
    _h5py as h5py
import numpy as np
import os
import shutil
import tempfile

from yt.analysis_modules.cosmological_observation.api import \
     LightCone
from yt.testing import \
    assert_equal, \
    requires_module
from yt.utilities.answer_testing.framework import \
    AnswerTestingTest, \
    requires_sim

ETC = "enzo_tiny_cosmology/32Mpc_32.enzo"

@requires_module("h5py")
class LightConeProjectionTest(AnswerTestingTest):
    _type_name = "LightConeProjection"
    _attrs = ()
 
    def __init__(self, parameter_file, simulation_type):
        self.parameter_file = parameter_file
        self.simulation_type = simulation_type
        self.ds = os.path.basename(self.parameter_file)

    @property
    def storage_name(self):
        return os.path.basename(self.parameter_file)

    def run(self):
        # Set up in a temp dir
        tmpdir = tempfile.mkdtemp()
        curdir = os.getcwd()
        os.chdir(tmpdir)

        lc = LightCone(
            self.parameter_file, self.simulation_type, 0., 0.1,
            observer_redshift=0.0, time_data=False)
        lc.calculate_light_cone_solution(
            seed=123456789, filename="LC/solution.txt")
        lc.project_light_cone(
            (600.0, "arcmin"), (60.0, "arcsec"), "density",
            weight_field=None, save_stack=True)

        fh = h5py.File("LC/LightCone.h5")
        data = fh["density_None"].value
        units = fh["density_None"].attrs["units"]
        assert units == "g/cm**2"
        fh.close()

        # clean up
        os.chdir(curdir)
        shutil.rmtree(tmpdir)

        mean = data.mean()
        mi = data[data.nonzero()].min()
        ma = data.max()
        return np.array([mean, mi, ma])

    def compare(self, new_result, old_result):
        assert_equal(new_result, old_result, verbose=True)

@requires_sim(ETC, "Enzo")
def test_light_cone_projection():
    yield LightConeProjectionTest(ETC, "Enzo")
