import os
import sys
import numpy as np
import igraph as ig
import pandas as pd
import copy

def personalized_pagerank(P, i, alpha, epsilon=1e-9, warnings=False):
    """
    This function is used to determine the vector probability using a PPR approach applied on the graph without the targeted node, only considering its neighbours.
    In this second approach we do not remove the targeted node to prevent the creation of pseudo sub-components while the damping factor get close to 1.
    We compute the PPR using all the network, also restarting from the neighborhood of the target compound.
    Then, we re-estimate the probability vector, setting the probability to be on the targeted node at 0, to represent the proportion of time passed out of the targeted node.
    Args:
        P ([numpy.ndarray]): Transition probability matrix of the graph
        i ([int]): Index of the target node
        alpha (float): The damping factor. WARGNING alpha is [0, 1[. '1' is excluded because we need restart probabilies to ensure the graph connexion !
        epsilon ([float], optional): Tolerance for convergence. Defaults to 1e-9.

    Returns:
        [numpy.ndarray]: Vector of probabilities to be out of the targeted node
    """
    # Get length
    l = P.shape[0]

    # Create restart vector on the targeted node
    v = np.zeros((l, 1))
    v[i, 0] = 1

    # No need to deal with sink nodes as we only consider connected graphs
    e = np.ones((l, 1))
    M = alpha * P + ((1 - alpha) * e) @ v.T

    # Apply Power method
    c = 1
    M = M.T
    pi = v
    new_pi = M @ v
    while np.linalg.norm(pi - new_pi) > epsilon:
        pi = new_pi
        new_pi = M @ pi
        c += 1

    # Float are basically imprecise and so after several matrix multiplications, the sum of probabilities in the vector may not equal to 1, but 0.9999999999999999 or 1.0000000000000001 for example. 
    if np.sum(new_pi, axis=0, dtype=np.float16) != 1 and warnings:
        print("[WARNING] At index " + str(i) + ": the final probability vector does not sum to 1. This may be due to float approximation errors")

    return new_pi

def transition_probability_matrix(g, alpha):
    """This function is used to compute the PPR for all nodes.

    Args:
        g (igraph.Graph): The compound graph
        alpha (float, optional): The damping factor.

    Returns:
        numpy array: The probability matrix. Probability vectors are stored in columns ! For instance, P_{A,B} = Probability that a mention starting from B reach A
    """

    # Get adjacency matrix
    if g.is_directed():

        # Use weight of edges as probabilities
        P = np.array(g.get_adjacency(attribute='Weight').data)

        # Test that all proba sum to 1. Use np.round(x, 9) to avoid float approximation errors
        if not sum(np.round(P.sum(axis=1), 9) == 1) == g.vcount():
            print("[ERROR] All the probability in the transition matrix extracted from edges' weights don't sum to 1 !")
            sys.exit(1)

    else:
        A = np.array(g.get_adjacency().data)
        P = np.diag(1/A.sum(axis=1)) @ A

    # If alpha is set to 0, simply return the probability matrix from direct neighborhood in columns using transpose, otherwise compute PPR
    if not alpha:
            full = P.T
    else:
        full = np.zeros(P.shape)
        for i in range(0, P.shape[0]):
            full[:, i] = personalized_pagerank(P, i, alpha)[:, 0]

    return full


def compute_weights(probabilities, table_species_corpora, q):
    """This function is used to compute the weight matrix. The weight matrix contains in columns the contributions of each compounds on the prior on the targeted compound. 
    For instance, W_{A,B} = proportion (probability) that a mention arrived on B, comes from A

    Args:
        probabilities (np.array): The PPR pobability matrix (Cf. transition_probability_matrix)
        table_species_corpora (pandas.DataFrame): A Dataframe containing the corpus size of each compounds
        q (float): The tolerance threshold. If the probability P_{A,B}, that a mention starting from B reach A is lower than the tolerance threshold, we consider that the compound B is not legitimate to contribute to the prior of A.

    Returns:
        [np.array]: The weight matrix
    """
    sigmas = copy.copy(probabilities)

    # To determine contributors of each compounds, we build a constrain matrix using the q parameter (tolerance threshold) to filter out contributors that a too far from the targeted node and also to avoid a compound to contribute itself to its prior.
    # First we get probability from the PPR
    constrains = copy.copy(probabilities)

    # To filter using the tolerance threshold, we first re-estimate probability without considering self-contribution !
    np.fill_diagonal(constrains, 0)
    constrains = constrains @ np.diag(1/(constrains.sum(axis=0)))

    # Test which resulting probability are higher than the tolerance threshold to create the constrain matrix. The diagonal will always be at False as we set 0 on diag previously
    constrains = (constrains > q) * 1

    # Apply constrain matrix on probability
    sigmas = sigmas * constrains

    # We re-estimate probabilities after filtering from constrains
    sigmas = sigmas @ np.diag(1/(sigmas.sum(axis=0)))

    # Get vector of corpora sizes
    v = np.array([table_species_corpora["TOTAL_PMID_SPECIE"]]).T

    # Compute totals
    t = sigmas @ v

    # But sometimes, the neighborhood does not have any mentions to transmit and the total recieved may be 0. to avoid divide by 0 we add 1 to every time there is 0.
    t = t + (t == 0) * 1

    # Compute value by specie
    w = (sigmas @ np.diag(v[:, 0])).T

    # Normalise by total
    weights = w @ np.diag(1/t[:, 0])

    # Weights are also store in columns !
    return weights


def create_probabilities_and_weights(g, g_name, alpha, table_species_corpora, q, id_att):
    """
    This function is used to achieve the creation of the probability and weight tables, by computing the PPR on the compound graph, applying the threshold for the neighborhood of influence and then calculating the associated weights.
    Also, this function create a cache directory to store probability and weight tables for each already computed alpha, avoiding to re-compute them.

    Args:
        g (igraph.Graph): The compound graph
        g_name (str): Name of the graph used to create the cache directory with probabilities and weights. By default the file name will be used.
        alpha (float): The damping factor.
        table_species_corpora (pandas.DataFrame): table of specie corpora
        q (float): The tolerance threshold for neighborhood influence
        id_att (str): The name of the vertex attribute containing the specie identifier (eg. M_m0001c) in the graph.

    Returns:
        [np.array, np.array]: array for probabilities and weights
    """

    cache_proba_dir_path = os.path.join("./cache", g_name,"PROBA")
    cache_weights_dir_path = os.path.join("./cache", g_name,"WEIGHTS")

    # If cache dir does not exists, create and fill it:
    for path in [cache_proba_dir_path, cache_weights_dir_path]:

        if not os.path.exists(path):
            try:
                os.makedirs(path, 0o755, exist_ok=False)

            except OSError:
                print("[INFO] Creation of the directory %s failed" % path)

            else:
                print("[INFO] Successfully created the directory %s" % path)

    # Test if probabilities for required alpha as already been computed :
    proba_path = os.path.join(cache_proba_dir_path, "PROBA_" + str(alpha) + ".csv")

    if not os.path.exists(proba_path):

        # Compute probabilities
        print("\n[INFO] Compute probabilities using alpha = " + str(alpha))
        probabilities = transition_probability_matrix(g, alpha=alpha)

        # Round probabilities with 9 decimals:
        probabilities = np.around(probabilities, 9)

        # Write probabilities in cache:
        df_probabilities = pd.DataFrame(probabilities, columns=g.vs[id_att], index=g.vs[id_att])
        df_probabilities.to_csv(proba_path, index=True, header=True)

    else:
        print("\n[INFO] Get probabilities with alpha = " + str(alpha) + " in cache dir")
        probabilities = pd.read_csv(proba_path, index_col=0, header=0).to_numpy()

    # Test if probabilities for required alpha as already been computed :
    weight_path = os.path.join(cache_weights_dir_path, "WEIGHTS_" + str(alpha) + ".csv")

    if not os.path.exists(weight_path):
        print("\n[INFO] Compute weights using alpha = " + str(alpha))
        weights = compute_weights(probabilities, table_species_corpora, q)

        # Round probabilities with 9 decimals:
        weights = np.around(weights, 9)

        # Write probabilities in cache:
        df_weights = pd.DataFrame(weights, columns=g.vs[id_att], index=g.vs[id_att])
        df_weights.to_csv(weight_path, index=True, header=True)

    else:
        print("\n[INFO] Get weights with alpha = " + str(alpha) + " in cache dir")
        weights = pd.read_csv(weight_path, index_col=0, header=0).to_numpy()


    return probabilities, weights