import numpy as np
import unittest

from quiver import *

class InertiaTensorTestCase(unittest.TestCase):
    def setUp(self):
        self.system = System('../test/iitensor.out')
        masses = { '1': 3.0,
                   '3': 7.0,
                   '6': 11.0 }
        self.systemIso = Isotopologue(self.system, masses)

    def tearDown(self):
        pass

    def test_iitensor(self):
        pass

if __name__ == '__main__':
    unittest.main()
