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
    coocurences = pd.merge(label_to_index, data, on = "SPECIE", how = "left")
    
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
            - FOT: FinishOnTarget, for each node, the resulting vector contains the probability that a walker on the targeted node comes from a particular node. The result of the StartFromTarget propagation in used to compute the FinishOnTarget probabilities.
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


def observation_uninformative_prior(k, n):
    """This function is used to build the posterior distribution using just an uninformative prior Beta(1,1)

    Args:
        k (integer): The coocurence for the MeSH 'M' for the targeted specie (number of success)
        n (integer): The corpus size for the targeted specie (number of trials)
    
    Returns:
    [collection]: A collection with:
        - alpha (int): alpha parameter related to the posterior Beta distribution: k + 1
        - beta (int): beta parameters related to the posterior Beta distribution: (n - k) + 1
        - x (list): probabilities
        - y (list): density 
        - mu (float): mean of the posterior distribution
    """
    # Posterior using unformative prior:
    r = collections.namedtuple("uninformativeprior", ["alpha", "beta", "x", "f", "mu"])
    x = np.arange(0, 1, 0.001).tolist()
    # Get distribution using uniformative prior: Beta(1,1)
    alpha = k + 1
    beta = (n - k) + 1
    y = ss.beta.pdf(x, a = alpha, b = beta)

    # Get mean of the distribution (Bayes estimator of p)
    mu = (k + 1)/(n + 2)

    res = r(alpha, beta, x, y, mu)

    return res

def create_prior_beta_mix(weights, cooc , corpora):
    """This function is used to determine values of the prior mixture distribution.
    In the prior mixture distribution, individual components are Beta() distributions related to the probability 'p' of success: an article discussing about a specie 's', also discusses the MeSH descriptor 'M'  
    Weights used in the mixture model are probabilities that a walker on the targeted specie comes from a particular specie 's': FinishOnTarget

    Args:
        weights (list): A list of float used as weights in the mixture distribution. (Cf. propagation_volume mode FinishOnTarget)
        cooc (list): A list of integer values representing co-occurences between species in the graph and a particular MeSH descriptor  
        corpora (list):  A list of integer values representing copus sizes related to each compounds in the graph

    Returns:
        [collection]: A collection with:
            - alpha (list): vector of alpha parameters related to each individual Beta distribution in the prior mixture. 
            - beta (list): vector the beta parameters related to each individual Beta distribution in the prior mixture.
            - weights (list): vector of weights related to each individual Beta distribution in the prior mixture.
            - x (list): probabilities
            - y (list): density 
    """
    # Get parameters
    r = collections.namedtuple("priormix", ["alpha", "beta", "weights", "x", "f", "mu"])
    x = np.arange(0, 1, 0.001).tolist()
    l = len(weights)

    # Get beta component paremeters for each compounds, using a posterior with uniformative prior
    alpha = [(cooc[it] + 1) for it in range(0, l)]
    beta = [(corpora[it] - cooc[it] + 1) for it in range(0, l)]

    f_i = [ss.beta.pdf(x, a = alpha[it], b = beta[it]) for it in range(0, l)]
    # Create Beta mix:
    y = np.dot(weights, f_i)

    # Get mean
    mu_i = [(alpha[it]/(alpha[it] + beta[it])) for it in range(0, l)]
    mu = np.dot(weights, mu_i)

    mix = r(alpha, beta, weights, x, y, mu)

    return mix
    


def create_posterior_beta_mix(k, n, weights_pior, alpha_prior, beta_prior, use_log = True):
    """This function is used to compute the posterior mixture distribution. The prior mixture model is updated for the observation of the coocurences on the targeted specie

    Args:
        k (integer): The coocurence for the MeSH 'M' for the targeted specie (number of success)
        n (integer): The corpus size for the targeted specie (number of trials)
        weights_pior (list): The prior mixture distribution weights (item 'weights' in create_prior_beta_mix result)
        alpha_prior (list): The prior mixture distribution alpha parameters (item 'alpha' in create_prior_beta_mix result)
        beta_prior (list): The prior mixture distribution beta parameters (item 'beta' in create_prior_beta_mix result)
        use_log (boolean, optional): A boolean telling if the computation of the weights W have to be achieve using classic formula or using logs (As alpha and beta values are often large, log method is prefered). Default = True

    Returns:
        [collection]: A collection with:
            - alpha (list): vector of alpha parameters related to each individual Beta distribution in the posterior mixture. 
            - beta (list): vector the beta parameters related to each individual Beta distribution in the posterior mixture.
            - weights (list): vector of weights related to each individual Beta distribution in the posterior mixture.
            - x (list): probabilities
            - y (list): density 
    """
    r = collections.namedtuple("posteriormix", ["alpha", "beta", "weights", "x", "f", "mu"])
    l = len(weights_pior)
    x = np.arange(0, 1, 0.001).tolist()

    # Get posterior parameters
    alpha_post = [(alpha_prior[it] + k) for it in range(0, l)]
    beta_post = [(beta_prior[it] + (n - k)) for it in range(0, l)]

    # Get Weight
    if use_log:
        # Compute log of W_i. Indeed sc.beta goes to 0 when alpha and beta are large which lead to a 0 division. use Log(Beta(a,b)) allow to compute W in thoose cases
        C = [sc.betaln(alpha_post[it], beta_post[it]) - sc.betaln(alpha_prior[it], beta_prior[it]) for it in range(0, l)]
        Z = [np.log(weights_pior[it]) + C[it] for it in range(0, l)]
        W = [np.exp((Z[it] - (sc.logsumexp(Z)))) for it in range(0, l)]
    else:    
        C = [sc.beta(alpha_post[it], beta_post[it])/sc.beta(alpha_prior[it], beta_prior[it]) for it in range(0, l)]
        W = [(weights_pior[it] * C[it]/(np.dot(weights_pior, C))) for it in range(0, l)]

    # Create posterior distribution by componennts
    f_post_i = [ss.beta.pdf(x, a = alpha_post[it], b = beta_post[it]) for it in range(0, l)]

    # Get posteriors probabilities
    y = np.dot(W, f_post_i)

    # Get mean
    mu_i = [(alpha_post[it]/(alpha_post[it] + beta_post[it])) for it in range(0, l)]
    mu = np.dot(W, mu_i)

    mix = r(alpha_post, beta_post, W, x, y, mu)
    
    return mix



####################
### Computations ###
####################

def plot_distributions(uninformative, prior_mix, posterior_mix):
    plt.plot(uninformative.x, uninformative.f, label = "Posterior with uninformative prior")
    plt.plot(prior_mix.x, prior_mix.f, label = "prior mix")
    plt.plot(posterior_mix.x, posterior_mix.f, label = "Posterior with prior mix")
    plt.title('Differences between Uninformative, mix prior and posterior mix distribution ')
    plt.legend()
    plt.show()

def computation(specie, mesh, table_cooc, table_corpora, matrix_proba):
    # Get cooc vector. It only contains species that have at least one article, need to left join.
    cooc = table_cooc[table_cooc["MESH"] == mesh][["index", "COOC"]]

    # Create table to resume needed data
    data = pd.merge(table_corpora, cooc, on = "index", how = "left").fillna(0)
    data["weights"] = matrix_proba[specie].tolist()
    print(data)

    # Remove values for targeted index
    index = int(data[data["SPECIE"] == specie]["index"])
    weights = data["weights"].tolist()
    del weights[index]
    cooc = data["COOC"].tolist()
    k = cooc.pop(index)
    corpora = data["TOTAL_PMID_SPECIE"].tolist()
    n = corpora.pop(index)

    # Uninformative
    obs = observation_uninformative_prior(k, n)
    # Prior mix:
    prior_mix = create_prior_beta_mix(weights, cooc, corpora)

    # Posterior mix:
    posterior_mix = create_posterior_beta_mix(k, n, prior_mix.weights, prior_mix.alpha, prior_mix.beta)

    print("> Mean uninformative: " + str(obs.mu))
    print("> Mean prior mix: " + str(prior_mix.mu))
    print("> Mean posterior_mix: " + str(posterior_mix.mu))

    plot_distributions(obs, prior_mix, posterior_mix)
    
    return data