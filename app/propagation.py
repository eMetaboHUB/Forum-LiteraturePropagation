import sys, os
import igraph as ig
import cairo
import pandas as pd
import numpy as np



def import_metabolic_network(path, undirected = True, format = "gml"):
    """This function is used to import an metabolic network, by default as undirected

    Args:
        path (str): path to the metabolic
        undirected (bool, optional): Graph directed or undirected. Defaults to True.
        format (str, optional): Graph format. Defaults to gml.
    """
    if not os.path.isfile(path):
        print("Error: Can't find a file at " + path)
        return None
    try:
        g = ig.read(path, format = format)
    except Exception as e:
        print("Error while reading graph file at " + path)
        print(e)
        return None
    print("Metabolic network has been imported successfully.")
    print("> Number of vertices: " + str(g.vcount()))
    print("> Number of edges: " + str(g.ecount()))
    if undirected:
        g.to_undirected()
    return g
