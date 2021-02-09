import sys, os
import igraph as ig
import cairo
import pandas as pd
import numpy as np
import copy
import collections
np.set_printoptions(suppress=True)



###################
###### Utils ######
###################

def import_metabolic_network(path, undirected = True, format = "gml"):
    """This function is used to import an metabolic network, by default as undirected

    Args:
        path (str): path to the metabolic
        undirected (bool, optional): Graph directed or undirected. Defaults to True.
        format (str, optional): Graph format. Defaults to gml.

    Returns:
        [igraph.Graph]: The compound graph
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

def import_corpora_sizes(path, g, name_att = "label"):
    """This function is used to import specie corpora sizes. The file containing corpora sizes must contains two columns: SPECIE (specie label) and TOTAL_PMID_SPECIE (corpus size)

    Args:
        path (tring): path to the species corpus size file
        g (igraph.Graph): the compound graph, imported using import_metabolic_network
        name_att (str, optional):  The name of the vertex attribute containing names. Defaults to "label".

    Returns:
        [pandas DataFrame]: a DataFrame containg: index of species in the compound grand (col index), specie label (SPECIE) and specie corpus size (TOTAL_PMID_SPECIE)
    """
    # Get label to index from graph 
    label_to_index = pd.DataFrame({"index": range(0, len(g.vs)), "SPECIE": g.vs[name_att]})
    # Read data
    if not os.path.isfile(path):
        print("Error: Can't find a file at " + path)
        return None
    try:
        data = pd.read_csv(path)
    except Exception as e:
        print("Error while reading graph file at " + path)
        print(e)
        return None
    if "SPECIE" not in data:
        print("Error: missing column 'specie' in " + path)
        return None
    if "TOTAL_PMID_SPECIE" not in data:
        print("Error: missing column 'TOTAL_PMID_SPECIE' in " + path)
        return None
    corpora_sizes = pd.merge(label_to_index, data)
    return corpora_sizes



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

def propagation_volume(g, name_att = "label", direction = "both"):
    """This function is used to compute the PPR, excluding the targeted node itself, for each node of the graph

    Args:
        g (igraph.Graph): The compound graph
        name_att (str, optional): The name of the vertex attribute containing names. Defaults to "label".
        direction (str, optional): The direction og random walks that will be used to compute probabilities:
            - SFT: StartFromTarget, for each node the resulting vector contains the probabilities to be on a particular compound node during the random walk starting from the targeted node.
            - FOT: FinishOnTarget, for each node, the resulting vector contains the probabilites that a walker on the targeted node comes from a particular node. The result of the forward propagation in used to compute the backward probabilities.
            - both: The both matrix probabilities are computed are returned

    Returns:
        collections.namedtuple: A named.tuple containing pandas DataFrame representing the probability matrix using SFT and/or FOT propagation.
    """
    # Init tuple
    r = collections.namedtuple("propagation", ["SFT", "FOT"])

    # Compute for each node SFT propagation
    A = np.array(g.get_adjacency().data)
    full = np.zeros(A.shape)
    for i in range(0, A.shape[0]):
        full[:, i] = compute_PR(A, i)[:, 0]
    
    # If SFT direction
    if direction == "SFT":
        df_SFT = pd.DataFrame(full, columns=g.vs[name_att], index=g.vs[name_att])
        result = r(df_SFT, None)

    # If backward direction
    if direction == "FOT":
        d = np.diag(1/full.sum(axis = 1))
        bkw = (full.T) @ d
        df_FOT = pd.DataFrame(bkw, columns=g.vs[name_att], index=g.vs[name_att])
        
    
    # If both
    if direction == "both":
        d = np.diag(1/full.sum(axis = 1))
        bkw = (full.T) @ d
        df_SFT = pd.DataFrame(full, columns=g.vs[name_att], index=g.vs[name_att])
        df_FOT = pd.DataFrame(bkw, columns=g.vs[name_att], index=g.vs[name_att])
        result = r(df_SFT, df_FOT)
    
    return result




