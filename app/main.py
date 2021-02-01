import argparse
from propagation import import_metabolic_network

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument("--graph", help="path to metabolic network compound graph", type = str, required = True, dest = 'g_path')


args = parser.parse_args()
g = import_metabolic_network(args.g_path)