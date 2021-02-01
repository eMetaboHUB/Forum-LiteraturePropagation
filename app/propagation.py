import sys, os
import igraph as ig
import cairo
import pandas as pd
import numpy as np
import copy



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



###################
### Propagation ###
###################


def compute_PR(A, i, alpha = 0.8):
    # Get truncated length
    l = A.shape[0]
    # Create restart vector by extracting probability, fromated as a column vector
    v = np.array([(A[i, :]/A[i, :].sum())])
    # Delete row and column associated with the targeted index
    v = np.delete(v, i, 1)
    truncate_A = np.delete(np.delete(A, i, 0), i, 1)
    # Sink node vector
    a = np.array([((truncate_A.sum(axis = 1) == 0) * 1)])
    # For sink nodes, the diagonal element is 0 and 1 for non-sink nodes. Adding the vector a to diagonal elements, ensure that we will not divide by 0 for sink nodes
    z = truncate_A.sum(axis = 1) + a
    d = np.diag(1/z[0])
    # Get probability matrix
    P = d @ truncate_A
    e = np.ones((1, l - 1))
    print((alpha * a + (1 - alpha) * e).T)
    P_final = alpha * P + (alpha * a + (1 - alpha) * e).T @ v
    return P_final

def propagation_volume(g):
    A = np.array(g.get_adjacency().data)
    return A