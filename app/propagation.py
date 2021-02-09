import sys, os
import igraph as ig
import cairo
import pandas as pd
import numpy as np
import copy
import collections
import scipy.special as sc
import scipy.stats as ss
import matplotlib.pyplot as plt
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

def import_table(path):
    """This function is a method to read a table. It check if the file exists, and handle exceptions while importing it with pandas.

    Args:
        path (string): path to the table

    Returns:
        (pandas.DataFrame): The imported table
    """
    # Test if file exists
    if not os.path.isfile(path):
        print("Error: Can't find a file at " + path)
        return None
    # Read table
    try:
        data = pd.read_csv(path)
    except Exception as e:
        print("Error while reading graph file at " + path)
        print(e)
        return None
    
    return data

def import_and_map_indexes(path, g, name_att = "label"):
    """This function is used to import tables related to species. It also add a new column containing the species index in the graph g.
    Args:
        path (string): path to the table to import
        g (igraph.Graph): The compound graph
        name_att (str, optional): The name of the vertex attribute containing names in the graph. Defaults to "label".

    Returns:
        (pandas.DataFrame): The imported table with a new column, index, indicating the species index in the graph
    """
    # Create a table to map species labels (SPECIE column) to indexes in the graph
    label_to_index = pd.DataFrame({"index": range(0, len(g.vs)), "SPECIE": g.vs[name_att]})
    data = import_table(path)

    # Merge
    coocurences = pd.merge(label_to_index, data)
    
    return coocurences

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


####################
## Beta functions ##
####################

def create_prior_beta_mix(p, cooc , corpora):
    # Get parameters
    r = collections.namedtuple("priormix", ["p", "f"])
    l = len(p)
    x = np.linspace(0, 1, 1000)

    # Get beta component paremeters for each compounds, using a posterior with uniformative prior
    alpha = [(cooc[it] + 1) for it in range(0, l)]
    beta = [(corpora[it] - cooc[it] + 1) for it in range(0, l)]
    print(alpha)
    print(beta)
    f_i = [ss.beta.pdf(x, a = alpha[it], b = beta[it]) for it in range(0, l)]
    # Create Beta mix:
    y = np.dot(p, f_i)
    mix = r(x,y)
    return mix
    



def create_posterior_beta_mix(k, n, p, cooc , corpora, l):