import argparse
from propagation import *

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument("--graph", help="path to metabolic network compound graph", type = str, required = True, dest = 'g_path')
parser.add_argument("--specie.corpora", help="path to the species corpus size file ", type = str, required = True, dest = 'specie_corpora_path')
parser.add_argument("--specie.cooc", help="path to the species MeSH co-occurences file ", type = str, required = True, dest = 'specie_mesh_path')
parser.add_argument("--mesh.corpora", help="path to the MeSH corpus size file ", type = str, required = True, dest = 'mesh_corpora_path')

args = parser.parse_args()
g = import_metabolic_network(args.g_path)

print("> Import species corpora sizes ... ", end = '')
species_corpora = import_and_map_indexes(args.specie_corpora_path, g)
print("Ok")

print("> Import species-MeSH co-occurences ... ", end = '')
coocurences = import_and_map_indexes(args.specie_mesh_path, g)
print("Ok")

print("> Import MeSH corpora sizes ... ", end = '')
mesh_corpora = import_table(args.mesh_corpora_path)
print("Ok")

print(species_corpora)
print(coocurences)
print(mesh_corpora)
# test = propagation_volume(g)
# test.SFT.to_csv("SFT.csv")
# test.FOT.to_csv("FOT.csv")


