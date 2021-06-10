import unittest
import pandas as pd
from app.propagation import *

class TestPropagationMethods(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.g = import_metabolic_network("tests/data/test_urea.gml")
        cls.table_specie = table_species_corpora = import_and_map_indexes("tests/data/species_cid_pmid_full.csv", cls.g)
        cls.probabilities = propagation_volume(cls.g, 0.8)
        cls.q = 1/(len(cls.g.vs) - 1)
        cls.weights = compute_weights(cls.probabilities, cls.table_specie, cls.q)
        cls.prior_mix_uninf = create_prior_beta_mix(weights = [0.5, 0.5], cooc = [1, 3], corpora = [4, 4], seq = 0.001, alpha_prior = 1, beta_prior = 1)
        cls.prior_mix = create_prior_beta_mix(weights = [0.5, 0.5], cooc = [1, 3], corpora = [4, 4], seq = 0.001, alpha_prior = 4, beta_prior = 10)
        cls.posterior_mix = create_posterior_beta_mix(k = 8, n = 10, weights_pior = cls.prior_mix_uninf.weights, alpha_prior = cls.prior_mix_uninf.alpha, beta_prior = cls.prior_mix_uninf.beta, seq = 0.001)
        # Random examples for simple priors
        cls.simple_prior = simple_prior(3, 4, seq = 0.001)
        cls.simple_posterior = simple_posterior(2, 3, 3, 4, seq = 0.001)

    def test_import(self):
        self.assertEqual(type(self.g), ig.Graph)

    def test_propagation_volume(self):
        SFT_ref = pd.read_csv("tests/data/SFT_ref.csv", index_col = 0)
        test = pd.DataFrame(self.probabilities, columns=self.g.vs["label"], index=self.g.vs["label"])
        pd.testing.assert_frame_equal(test, SFT_ref)
    
    def test_comput_weights(self):
        W_ref = pd.read_csv("tests/data/W_ref_0.8.csv", index_col = 0)
        test = pd.DataFrame(self.weights, columns=self.g.vs["label"], index=self.g.vs["label"])
        pd.testing.assert_frame_equal(test, W_ref)
    
    def test_prior_mix(self):
        self.assertEqual(round(self.prior_mix_uninf.f[775], 5), 1.13562)

    def test_prior_glm(self):
        self.assertEqual(round(self.prior_mix.f[775], 7), 0.0049983)
    
    def test_posterior(self):
        self.assertEqual(round(self.posterior_mix.f[775], 5), 3.45663)

    def test_simple_prior(self):
        self.assertEqual(round(self.simple_prior.f[775], 5), 0.41049)
    
    def test_simple_posterior(self):
        self.assertEqual(round(self.simple_posterior.f[775], 5), 0.58248)


if __name__ == '__main__':
    unittest.main()