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

    def test_propagation_volume(self):
        df_ref = pd.read_csv("tests/data/proba_ref.csv", index_col = 0)
        pd.testing.assert_frame_equal(self.propagation_volume, df_ref)


if __name__ == '__main__':
    unittest.main()