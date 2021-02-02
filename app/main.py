import argparse
from propagation import *

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument("--graph", help="path to metabolic network compound graph", type = str, required = True, dest = 'g_path')


args = parser.parse_args()
g = import_metabolic_network(args.g_path)
test = propagation_volume(g)
test.to_csv("tests/data/proba.csv")


