import unittest
import pandas as pd
from app.propagation import *

class TestPropagationMethods(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.g = import_metabolic_network("tests/data/test_urea.gml")
        cls.propagation_volume = propagation_volume(cls.g)
        cls.prior_mix_uninf = create_prior_beta_mix(weights = [0.5, 0.5], cooc = [1, 3], corpora = [4, 4], seq = 0.001)
        cls.prior_mix_glm = create_prior_beta_mix(weights = [0.5, 0.5], cooc = [1, 3], corpora = [4, 4], seq = 0.001, alpha_prior = 4, beta_prior = 10)
        cls.posterior_mix = create_posterior_beta_mix(k = 8, n = 10, weights_pior = cls.prior_mix_uninf.weights, alpha_prior = cls.prior_mix_uninf.alpha, beta_prior = cls.prior_mix_uninf.beta, seq = 0.001)

    def test_import(self):
        self.assertEqual(type(self.g), ig.Graph)

    def test_propagation_volume_SFT(self):
        SFT_ref = pd.read_csv("tests/data/SFT_ref.csv", index_col = 0)
        pd.testing.assert_frame_equal(self.propagation_volume.SFT, SFT_ref)

    def test_propagation_volume_FOT(self):
        FOT_ref = pd.read_csv("tests/data/FOT_ref.csv", index_col = 0)
        pd.testing.assert_frame_equal(self.propagation_volume.FOT, FOT_ref)
    
    def test_prior_mix(self):
        self.assertEqual(round(self.prior_mix_uninf.f[775], 5), 1.13562)

    def test_prior_glm(self):
        self.assertEqual(round(self.prior_mix_glm.f[775], 7), 0.0049983)
    
    def test_posterior(self):
        self.assertEqual(round(self.posterior_mix.f[775], 5), 3.45663)


if __name__ == '__main__':
    unittest.main()