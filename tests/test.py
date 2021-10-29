import unittest
import pandas as pd
import numpy as np
from app.propagation import *

class TestPropagationMethods(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):

        # The graph
        cls.g = import_metabolic_network("tests/data/Human1/1.7/init/Human-GEM_1.7_CarbonSkeletonGraph.gml", undirected = False)
        
        # Table species corpora
        cls.table_species_corpora = import_and_map("tests/data/Human1/1.7/init/species_pmids_Human1_1.7.csv", cls.g, "label")
        cls.table_species_corpora = cls.table_species_corpora.fillna(0)
        
        # Table species - MeSH
        cls.table_coocurences = import_and_map("tests/data/Human1/1.7/init/species_mesh_pmids_Human1_1.7.csv", cls.g, "label")
        cls.table_mesh_corpora = cls.table_coocurences.groupby('MESH', as_index=False)[['COOC']].sum().rename(columns={"COOC": "TOTAL_CPD_MENTION_MESH"})
        cls.table_mesh_corpora["P"] = cls.table_mesh_corpora["TOTAL_CPD_MENTION_MESH"]/(cls.table_species_corpora['TOTAL_PMID_SPECIE'].sum())

        # PROBA
        # cls.probabilities = propagation_volume(cls.g, 0.4)
        # cls.probabilities = np.around(cls.probabilities, 9)
        cls.probabilities = pd.read_csv("tests/data/Human1/1.7/validation/PROBA_0.4.csv", index_col = 0, header = 0).to_numpy()

        # WEIGHTS
        cls.q = 1/(len(cls.g.vs) - 1)
        # cls.weights = compute_weights(cls.probabilities, cls.table_species_corpora, cls.q)
        # cls.weights = np.around(cls.weights, 9)
        cls.weights = pd.read_csv("tests/data/Human1/1.7/validation/WEIGHTS_0.4.csv", index_col = 0, header = 0).to_numpy()

        # TEST PRIOR/POSTERIOR MIX FUNCTIONS
        
        cls.prior_mix = create_prior_beta_mix(weights = [0.5, 0.5], cooc = [1, 3], corpora = [4, 4], seq = 0.001, alpha_prior = 4, beta_prior = 10)
        cls.posterior_mix = create_posterior_beta_mix(k = 8, n = 10, weights_pior = cls.prior_mix.weights, alpha_prior = cls.prior_mix.alpha, beta_prior = cls.prior_mix.beta, seq = 0.001)
        # Random examples for simple priors
        cls.simple_prior = simple_prior(3, 4, seq = 0.001)
        cls.simple_posterior = simple_posterior(2, 3, 3, 4, seq = 0.001)

    def test_import(self):
        self.assertEqual(self.g.vcount(), 2742)
        self.assertEqual(self.g.ecount(), 10480)
    
    def test_table_species_corpora(self):
        ref = pd.read_csv("tests/data/Human1/1.7/validation/table_species_corpora.csv", index_col = False)
        pd.testing.assert_frame_equal(self.table_species_corpora, ref)
    
    def test_table_mesh_corpora(self):
        ref = pd.read_csv("tests/data/Human1/1.7/validation/table_mesh_corpora.csv", index_col = False)
        pd.testing.assert_frame_equal(self.table_mesh_corpora, ref)

    def test_probabilities(self):
        ref = pd.read_csv("tests/data/Human1/1.7/validation/PROBA_0.4.csv", index_col = 0, header = 0).to_numpy()
        np.testing.assert_array_equal(ref, self.probabilities)
    
    def test_weights(self):
        ref = pd.read_csv("tests/data/Human1/1.7/validation/WEIGHTS_0.4.csv", index_col = 0, header = 0).to_numpy()
        np.testing.assert_array_equal(ref, self.weights)
    
    def test_computation(self):
        index = 2521
        data = pd.read_csv("tests/data/Human1/1.7/init/data.csv", index_col = False)
        p = 0.002214000288021922
        alpha_prior = 2.2140002880219223
        beta_prior = 997.7859997119781
        res = computation(index, data, p, alpha_prior, beta_prior, seq = 0.0001, report = False)
        self.assertEqual(np.round(res['CDF'], 9), 0.002127064)
        self.assertEqual(np.round(res['Log2FC'], 9), 3.216560006)
        self.assertEqual(np.round(res['priorCDF'], 9), 0.008455974)
        self.assertEqual(np.round(res['priorLog2FC'], 9), 4.054079044)

    def test_intial_prior(self):
        self.assertEqual(round(self.prior_mix.f[775], 7), 0.0049983)
    
    def test_posterior(self):
        self.assertEqual(round(self.posterior_mix.f[775], 5), 0.10279)

    def test_simple_prior(self):
        self.assertEqual(round(self.simple_prior.f[775], 5), 0.41049)
    
    def test_simple_posterior(self):
        self.assertEqual(round(self.simple_posterior.f[775], 5), 0.58248)
    
if __name__ == '__main__':
    unittest.main()