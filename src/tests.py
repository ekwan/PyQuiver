import numpy as np
import unittest

from quiver import *
'''
class InertiaTensorTestCase(unittest.TestCase):
    def setUp(self):
        self.system = System('../test/iitensor.out')
        masses = { 0: 3.0,
                   1: 7.0,
                   2: 11.0 }
        self.systemIso = Isotopologue(self.system, masses)

    def tearDown(self):
        pass

    def test_iitensor(self):
        pass
'''
class RotationTestCase(unittest.TestCase):
    def setUp(self):
        self.system = System('../test/rotation3.out')
        masses = { 0: 1.0,
                   1: 4.0,
                   2: 9.0,
                   3: 11.0,
                   4: 1.0,
                   5: 1.0, }
        self.systemIso = Isotopologue(self.system, masses)

    def tearDown(self):
        pass

    def test_iitensor(self):
        pass

if __name__ == '__main__':
    unittest.main()
