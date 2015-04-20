from __future__ import absolute_import
from yt.config import ytcfg
ytcfg["yt","loglevel"] = '50'
ytcfg["yt","suppressStreamLogging"] = 'True'

from . import hydro_tests # Just importing will register the tests!
from . import halo_tests
from . import particle_tests
from .runner import RegressionTestRunner

first_runner = RegressionTestRunner("first")
first_runner.run_all_tests()
second_runner = RegressionTestRunner("second", "first")
second_runner.run_all_tests()
