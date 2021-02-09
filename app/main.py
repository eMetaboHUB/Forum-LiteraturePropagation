import argparse
from propagation import *

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument("--graph", help="path to metabolic network compound graph", type = str, required = True, dest = 'g_path')
parser.add_argument("--specie.corpora", help="path to the species corpus size file ", type = str, required = True, dest = 'specie_corpora_path')

args = parser.parse_args()
g = import_metabolic_network(args.g_path)
corpora_sizes = import_corpora_sizes(args.specie_corpora_path, g)
print(corpora_sizes)
# test = propagation_volume(g)
# test.SFT.to_csv("SFT.csv")
# test.FOT.to_csv("FOT.csv")


