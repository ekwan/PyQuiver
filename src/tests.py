import numpy as np
import unittest

from quiver import *

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

if __name__ == '__main__':
    unittest.main()
