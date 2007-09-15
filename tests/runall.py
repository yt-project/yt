"""
Basic unit test suite to run all the tests at once. Hopefully this it's
clear how to append additional tests.

You can run this using: 

$ python runall.py

This can be done from any directory.

Currently, YT needs to be installed to run the tests, so that Python
can find the modules.
"""

import unittest

from test_lagos import TestLagos
from test_raven import TestRaven

def get_suite():
    suite_l = unittest.TestLoader().loadTestsFromTestCase(TestLagos)
    suite_r = unittest.TestLoader().loadTestsFromTestCase(TestRaven)
    suite = unittest.TestSuite([suite_l, suite_r])
    return suite

if __name__ == '__main__':
    suite = get_suite()
    unittest.main()


