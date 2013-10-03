"""
API for enzo_test



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from .runner import \
    RegressionTestRunner, \
    RegressionTestStorage, \
    clear_registry, \
    registry_entries

from .output_tests import \
    DatasetTest, \
    create_test

from .default_tests import \
    TestFieldStatistics, \
    TestAllProjections

from .xunit import \
    Xunit

from .halo_tests import \
    TestHaloCompositionHashHOP, \
    TestHaloCompositionHashFOF, \
    TestHaloCompositionHashPHOP

from .boolean_region_tests import \
    TestBooleanANDGridQuantity, \
    TestBooleanORGridQuantity, \
    TestBooleanNOTGridQuantity, \
    TestBooleanANDParticleQuantity, \
    TestBooleanORParticleQuantity, \
    TestBooleanNOTParticleQuantity

try:
    from .framework import AnswerTesting
except ImportError:
    raise

def run_nose(verbose=False, run_answer_tests=False, answer_big_data=False):
    import nose, os, sys, yt
    from yt.funcs import mylog
    orig_level = mylog.getEffectiveLevel()
    mylog.setLevel(50)
    nose_argv = sys.argv
    nose_argv += ['--exclude=answer_testing','--detailed-errors']
    if verbose:
        nose_argv.append('-v')
    if run_answer_tests:
        nose_argv.append('--with-answer-testing')
    if answer_big_data:
        nose_argv.append('--answer-big-data')
    initial_dir = os.getcwd()
    yt_file = os.path.abspath(yt.__file__)
    yt_dir = os.path.dirname(yt_file)
    os.chdir(yt_dir)
    try:
        nose.run(argv=nose_argv)
    finally:
        os.chdir(initial_dir)
        mylog.setLevel(orig_level)
