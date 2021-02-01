import igraph as ig
import cairo
import pandas as pd
import numpy as np
import sys



def import_metabolic_network(path, undirected = True, format = "gml"):
    """This function is used to import an metabolic network, by default as undirected

    Args:
        path (str): path to the metabolic
        undirected (bool, optional): Graph directed or undirected. Defaults to True.
        format (str, optional): Graph format. Defaults to gml.
    """
    try:
        g = ig.read(path, format = format)
    except FileNotFoundError as e:
        print("Unable to read source file at " + path)
        print(e)
    print("Metabolic network has been imported successfully.")
    print("> Number of vertices: " + str(g.vcount()))
    print("> Number of edges: " + str(g.ecount()))
    g.to_undirected()
    return g
