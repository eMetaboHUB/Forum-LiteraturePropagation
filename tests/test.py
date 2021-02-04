import unittest
import pandas as pd
from app.propagation import *

class TestPropagationMethods(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.g = import_metabolic_network("tests/data/test_urea.gml")
        cls.propagation_volume = propagation_volume(cls.g)

    def test_import(self):
        self.assertEqual(type(self.g), ig.Graph)

    def test_propagation_volume_SFT(self):
        SFT_ref = pd.read_csv("tests/data/SFT_ref.csv", index_col = 0)
        pd.testing.assert_frame_equal(self.propagation_volume.SFT, SFT_ref)

    def test_propagation_volume_FOT(self):
        FOT_ref = pd.read_csv("tests/data/FOT_ref.csv", index_col = 0)
        pd.testing.assert_frame_equal(self.propagation_volume.FOT, FOT_ref)


if __name__ == '__main__':
    unittest.main()