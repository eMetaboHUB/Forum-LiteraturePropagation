import argparse
from propagation import *

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument("--graph", help="path to metabolic network compound graph", type = str, required = True, dest = 'g_path')


args = parser.parse_args()
g = import_metabolic_network(args.g_path)
test = propagation_volume(g)
test.SFT.to_csv("SFT.csv")
test.FOT.to_csv("FOT.csv")


