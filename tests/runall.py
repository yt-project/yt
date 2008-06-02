"""
Basic unit test suite to run all the tests at once. Hopefully this it's
clear how to append additional tests.

You can run this using: 

$ python tests/runall.py

This should be done from the root directory of the installation.

YT can either be installed globally, or the extensions build with:

$ python setup.py build_ext --inplace
"""

import unittest

import test_lagos
import test_raven
import test_hdf5_reader

def get_suite():
    suite_l = unittest.defaultTestLoader.loadTestsFromModule(test_lagos)
    suite_r = unittest.defaultTestLoader.loadTestsFromModule(test_raven)
    suite_h = unittest.defaultTestLoader.loadTestsFromModule(test_hdf5_reader)
    suite = unittest.TestSuite([suite_l, suite_r])
    return suite

if __name__ == '__main__':
    unittest.main(defaultTest='get_suite')
