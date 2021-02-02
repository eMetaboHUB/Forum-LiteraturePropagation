import argparse
from propagation import *

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument("--graph", help="path to metabolic network compound graph", type = str, required = True, dest = 'g_path')


args = parser.parse_args()
g = import_metabolic_network(args.g_path)
A = propagation_volume(g)

c0 = compute_PR(A, 0)
c1 = compute_PR(A, 1)

print(c1.shape)
print(np.sum(c1, axis = 0))

print(c0.shape)
print(np.sum(c0, axis = 0))
