import argparse
from propagation import *

# Get arguments
# python app/main.py --graph="tests/data/test_urea.gml" --specie.corpora="data/species_cid_pmid.csv" --specie.cooc="data/species_cid_mesh_pmid.csv" --mesh.corpora="data/mesh_pmid.csv"
parser = argparse.ArgumentParser()
parser.add_argument("--graph", help="path to metabolic network compound graph", type = str, required = True, dest = 'g_path')
parser.add_argument("--specie.corpora", help="path to the species corpus size file ", type = str, required = True, dest = 'specie_corpora_path')
parser.add_argument("--specie.cooc", help="path to the species MeSH co-occurences file ", type = str, required = True, dest = 'specie_mesh_path')
parser.add_argument("--mesh.corpora", help="path to the MeSH corpus size file ", type = str, required = True, dest = 'mesh_corpora_path')

args = parser.parse_args()
g = import_metabolic_network(args.g_path)

print("> Import species corpora sizes ... ", end = '')
table_species_corpora = import_and_map_indexes(args.specie_corpora_path, g)
print("Ok")

print("> Import species-MeSH co-occurences ... ", end = '')
table_coocurences = import_and_map_indexes(args.specie_mesh_path, g)
print("Ok")

print("> Import MeSH corpora sizes ... ", end = '')
table_mesh_corpora = import_table(args.mesh_corpora_path)
print("Ok")


probabilities = propagation_volume(g)

mesh = "D056806"
specie = "M_urea"

table_test = computation(specie, mesh, table_coocurences, table_species_corpora, probabilities.FOT)



# plt.plot(prior_test.x, prior_test.f)
# plt.show()

