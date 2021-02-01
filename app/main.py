import argparse
from propagation import *

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument("--graph", help="path to metabolic network compound graph", type = str, required = True, dest = 'g_path')


args = parser.parse_args()
g = import_metabolic_network(args.g_path)
A = propagation_volume(g)
c = compute_PR(A, 1)
print(g.vs["label"])
print(c.shape)
print(c.sum(axis = 1))
print(c)