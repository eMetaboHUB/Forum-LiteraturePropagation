import sys, os
import igraph as ig
import cairo
import pandas as pd
import numpy as np
import copy
np.set_printoptions(suppress=True)


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
    print("> Metabolic network has been imported successfully.")
    print("> Number of vertices: " + str(g.vcount()))
    print("> Number of edges: " + str(g.ecount()))
    if undirected:
        g.to_undirected()
    return g



###################
### Propagation ###
###################

def compute_PR(A, i, alpha = 0.8, epsilon = 1e-9):
    """
    This function is used to determine the vector probability using a PPR approach applied on the graph without the targeted node, only considering its neighbours.
    Args:
        A ([numpy.ndarray]): Graph adjacency matrix
        i ([int]): Index of the target node
        alpha (float, optional): The damping factor. Defaults to 0.8.
        epsilon ([float], optional): Tolerance for convergence. Defaults to 1e-9.

    Returns:
        [numpy.ndarray]: Vector of stationary probabilities (as column vector)
    """
    # Get truncated length
    l = A.shape[0]
    # Create restart vector by extracting probability, fromated as a column vector.
    v = np.array([(A[i, :]/A[i, :].sum())]).T

    # Delete row and column associated with the targeted index
    v = np.delete(v, i, 0)
    truncate_A = np.delete(np.delete(A, i, 0), i, 1)
    
    # Sink node vector, as column vector
    a = np.array([((truncate_A.sum(axis = 1) == 0) * 1)]).T
    
    # For sink nodes, the diagonal element is 0 instead of 1 for non-sink nodes. Adding the vector a to diagonal elements, ensure that we will not divide by 0 for sink nodes
    z = truncate_A.sum(axis = 1) + a.T
    d = np.diag(1/z[0])
    
    # Get probability matrix
    P = d @ truncate_A
    e = np.ones((l - 1, 1))
    M = alpha * P + (alpha * a + (1 - alpha) * e) @ v.T
    
    # Apply Power method
    # Use transpose of M in power method
    c = 1
    M = M.T
    pi = v
    new_pi = M @ v
    while(np.linalg.norm(pi - new_pi) > epsilon):
        pi = new_pi
        new_pi = M @ pi
        c += 1
    print(str(c) + " iterations to convergence.")

    # Insert 0 at targeted index
    r = np.insert(new_pi, i, 0, axis = 0)
    # Float are basically imprecise and so after several matrix multiplications, the sum of probabilities in the vector may not equal to 1, but 0.9999999999999999 or 1.0000000000000001 for example. 
    if np.sum(r, axis = 0, dtype = np.float16) != 1:
        print("Warning at index " + str(i) + ": the final probability vector does not sum to 1. This may be due to float approximation errors")
    
    return r

def propagation_volume(g, name_att = "label"):
    """This function is used to compute the PPR, excluding the targeted node itself, for each node of the graph

    Args:
        g ([igraph.Graph]): The compound graph
        name_att (str, optional): The name of the vertex attribute containing names. Defaults to "label".

    Returns:
        [pandas.DataFrame]: A data.frame representing the probability matrix obtained after propagation. Propability vectors are stored in columns.
    """
    A = np.array(g.get_adjacency().data)
    full = np.zeros(A.shape)
    for i in range(0, A.shape[0]):
        full[:, i] = compute_PR(A, i)[:, 0]
    df = pd.DataFrame(full, columns=g.vs[name_att], index=g.vs[name_att])
    return df