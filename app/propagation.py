import sys, os
import igraph as ig
import cairo
import pandas as pd
import numpy as np
import copy
import collections
import scipy.special as sc
import scipy.stats as ss
import progressbar
import matplotlib.pyplot as plt
np.set_printoptions(suppress=True)



###################
###### Utils ######
###################

def import_metabolic_network(path, undirected = True, format = "gml", largest_comp = True):
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
        print("> Used as undirected")
        g.to_undirected()
    if largest_comp:
        print("> Extract largest component")
        g = g.clusters().giant()
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
    # Create a table to map species labels (SPECIE column) to indexes in the graph.
    label_to_index = pd.DataFrame({"index": range(0, len(g.vs)), "SPECIE": g.vs[name_att]})
    data = import_table(path)

    # Merge
    coocurences = pd.merge(label_to_index, data, on = "SPECIE", how = "left")
    
    return coocurences

###################
### Propagation ###
###################

def compute_PR(A, i, alpha, epsilon = 1e-9):
    """
    This function is used to determine the vector probability using a PPR approach applied on the graph without the targeted node, only considering its neighbours.
    Args:
        A ([numpy.ndarray]): Graph adjacency matrix
        i ([int]): Index of the target node
        alpha (float): The damping factor. WARGNING alpha is [0, 1[. '1' is excluded because we need restart probabilies to ensure the graph connexion !
        epsilon ([float], optional): Tolerance for convergence. Defaults to 1e-9.

    Returns:
        [numpy.ndarray]: Vector of stationary probabilities (as column vector)
    """
    # Get length
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
    # print(str(c) + " iterations to convergence.")

    # Insert 0 at targeted index
    r = np.insert(new_pi, i, 0, axis = 0)
    # Float are basically imprecise and so after several matrix multiplications, the sum of probabilities in the vector may not equal to 1, but 0.9999999999999999 or 1.0000000000000001 for example. 
    if np.sum(r, axis = 0, dtype = np.float16) != 1:
        print("Warning at index " + str(i) + ": the final probability vector does not sum to 1. This may be due to float approximation errors")

    return r

def compute_PR_2(A, i, alpha, epsilon = 1e-9):
    """
    This function is used to determine the vector probability using a PPR approach applied on the graph without the targeted node, only considering its neighbours.
    In this second approach we do not remove the targeted node to prevent the creation of pseudo sub-components while the damping factor get close to 1.
    We compute the PPR using all the network, also restarting from the neighborhood of the target compound.
    Then, we re-estimate the probability vector, setting the probability to be on the targeted node at 0, to represent the proportion of time passed out of the targeted node.
    Args:
        A ([numpy.ndarray]): Graph adjacency matrix
        i ([int]): Index of the target node
        alpha (float): The damping factor. WARGNING alpha is [0, 1[. '1' is excluded because we need restart probabilies to ensure the graph connexion !
        epsilon ([float], optional): Tolerance for convergence. Defaults to 1e-9.

    Returns:
        [numpy.ndarray]: Vector of probabilities to be out of the targeted node
    """
    #TODO NEXT FOR PROBA MATRIX, USE: g.get_adjacency(attribute='weight') pour extraire la matrice d'adj avec les poids (proba)
    # Get length
    l = A.shape[0]
    # Create restart vector on the targeted node
    v = np.zeros((l, 1))
    v[i,0] = 1

    # Sink node vector, as column vector
    a = np.array([((A.sum(axis = 1) == 0) * 1)]).T
    
    # For sink nodes, the diagonal element is 0 instead of 1 for non-sink nodes. Adding the vector a to diagonal elements, ensure that we will not divide by 0 for sink nodes
    # For the time, we use undirected graphs so there is no reason to have sink nodes, but it could happens if we use directed graphs
    z = A.sum(axis = 1) + a.T
    d = np.diag(1/z[0])
    
    # Get probability matrix
    P = d @ A
    e = np.ones((l, 1))
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
    # print(str(c) + " iterations to convergence.")

    # Float are basically imprecise and so after several matrix multiplications, the sum of probabilities in the vector may not equal to 1, but 0.9999999999999999 or 1.0000000000000001 for example. 
    if np.sum(new_pi, axis = 0, dtype = np.float16) != 1:
        print("Warning at index " + str(i) + ": the final probability vector does not sum to 1. This may be due to float approximation errors")

    return new_pi

def propagation_volume(g, alpha = 0.8, name_att = "label"):
    """This function is used to compute the PPR, excluding the targeted node itself, for each node of the graph

    Args:
        g (igraph.Graph): The compound graph
        alpha (float, optional): The damping factor. Defaults to 0.8.
        name_att (str, optional): The name of the vertex attribute containing names. Defaults to "label".

    Returns:
        numpy array: The probability matrix
    """
    # Get adjacency matrix
    A = np.array(g.get_adjacency().data)
    
    # If alpha is set to 0, simply return the probability matrix from direct neighborhood, otherwise compute PPR
    if not alpha:
        full = A @ np.diag(1/A.sum(axis=0))
    else:
        full = np.zeros(A.shape)
        for i in range(0, A.shape[0]):
            full[:, i] = compute_PR_2(A, i, alpha)[:, 0]

    df = pd.DataFrame(full, columns=g.vs[name_att], index=g.vs[name_att])
    
    return df


def compute_weights(probabilities, table_species_corpora, q):
    # Compute weights
    sigmas = probabilities.to_numpy()
    # To determine contributors of each compounds, we build a constrain matrix using the q parameter (tolerance threshold) to filter out contributors that a too far from the targeted node and also to avoid a compound to contribute itself to its prior.
    # First we get probability from the PPR
    constrains = copy.copy(sigmas)
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
    w = (sigmas @ np.diag(v[:,0])).T
    # Normalise by total
    weights = w @ np.diag(1/t[:, 0])
    
    # Weights are also store in columns !
    return weights

####################
## Beta functions ##
####################


def observation_uninformative_prior(k, n, seq, sampling = True):
    """This function is used to build the posterior distribution using just an uninformative prior Beta(1,1)

    Args:
        k (integer): The coocurence for the MeSH 'M' for the targeted specie (number of success)
        n (integer): The corpus size for the targeted specie (number of trials)
        seq (float): the step used in np.arange to create the x vector of probabilities.
    
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
    x = None
    y = None

    # Get distribution using uniformative prior: Beta(1,1)
    alpha = k + 1
    beta = (n - k) + 1
    if sampling:
        x = np.arange(0, 1 + seq, seq).tolist()
        y = ss.beta.pdf(x, a = alpha, b = beta)
    # Get mean of the distribution (Bayes estimator of p)
    mu = (k + 1)/(n + 2)

    res = r(alpha, beta, x, y, mu)

    return res

def estimate_prior_distribution_mesh_V2(mu, sample_size):
    """This function is used to estimate the prior distribution of the probability that a mention of a compound in an article, also involved the studied MeSH.
    This prior distribution is based on the independance hypothesis. 
    Args:
        mu (float): The expected mean of the prior probability distributions
        sample_size (integer): The expected sample size from which this prior distribution has been generated, related to the variability, certainity of the distribution. This parameter could also be refered as the concentration parameter. 

    Returns:
        [collection]: A collection with:
            - alpha: the alpha parameter of the prior distribution 
            - beta: the beta parameter of the prior distribution 
    """
    r = collections.namedtuple("prior_mesh", ["alpha", "beta"])

    alpha =  mu * sample_size
    beta = (1 - mu) * sample_size

    result = r(alpha, beta)
    return result

def estimate_prior_distribution_mesh(mesh_corpora):
    
    r = collections.namedtuple("prior_mesh", ["alpha", "beta"])
    
    if not mesh_corpora:
        result = r(1, 1)
        return result
    # Model parameters: 
    mu_intercept = -14.2634974
    mu_factor = 0.8363094
    sigma_intercept = -13.1256575
    sigma_factor = 0.7598162
    
    # Compute estimates:
    mu = sc.expit(mu_intercept + mu_factor * np.log(mesh_corpora))
    sigma = np.exp(sigma_intercept + sigma_factor * np.log(mesh_corpora))

    # Return parameters:
    alpha = mu/sigma
    beta = (1 - mu)/sigma

    result = r(alpha, beta)
    return result

def simple_prior(alpha_prior, beta_prior, seq, sampling = True):
    """This function is used to compute the density of a simple prior distribution, no mixture.

    Args:
        alpha_prior (float): The alpha parameter of the prior probability distribution (Cf. estimate_prior_distribution_mesh_V2)
        beta_prior (float): The beta parameter of the prior probability distribution (Cf. estimate_prior_distribution_mesh_V2)
        seq (float): The step used in np.arange to create the x vector of probabilities, only when sampling is True.
        sampling (bool, optional): Does the function have to compute a sampling of density values ?

    Returns:
    [collection]: A collection with:
        - alpha (float): The alpha parameter to the prior distribution. 
        - beta (float): The beta parameter to the prior distribution.
        - x (list): Probabilities
        - f (list): Densities
        - mu (float): Mean of the distribution  
    """
    # Get parameters
    r = collections.namedtuple("simple_prior", ["alpha", "beta", "x", "f", "mu"])
    x = None
    y = None
    
    if sampling:
        x = np.arange(0, 1 + seq, seq).tolist()
        y = ss.beta.pdf(x, a = alpha_prior, b = beta_prior)

    # Get mean
    mu = alpha_prior/(alpha_prior + beta_prior)
    res = r(alpha_prior, beta_prior, x, y, mu)

    return res

def simple_posterior(cooc, corpora, alpha_prior, beta_prior, seq, sampling = True):
    """This function is used to estimate parameters and density of a simple posterior distribution, no mixture.

    Args:
        cooc (list): The co-occurence between the specie in the graph and a particular MeSH descriptor  
        corpora (list): The corpus size related to the specie in the graph
        alpha_prior (float): The alpha parameter of the prior probability distribution (Cf. estimate_prior_distribution_mesh_V2)
        beta_prior (float): The beta parameter of the prior probability distribution (Cf. estimate_prior_distribution_mesh_V2)
        seq (float): The step used in np.arange to create the x vector of probabilities, only when sampling is True.
        sampling (bool, optional): Does the function have to compute a sampling of density values ?

    Returns:
    [collection]: A collection with:
        - alpha (float): The alpha parameter to the post distribution. 
        - beta (float): The beta parameter to the prior distribution.
        - x (list): Probabilities
        - f (list): Densities
        - mu (float): Mean of the distribution  
    """
    # Get parameters
    r = collections.namedtuple("simple_posterior", ["alpha", "beta", "x", "f", "mu"])
    x = None
    y = None

    # Get beta component paremeters for each compounds, using a posterior with uniformative prior
    alpha = cooc + alpha_prior
    beta = corpora - cooc + beta_prior
    
    if sampling:
        x = np.arange(0, 1 + seq, seq).tolist()
        y = ss.beta.pdf(x, a = alpha, b = beta)

    # Get mean
    mu = alpha/(alpha + beta)
    res = r(alpha, beta, x, y, mu)

    return res

def create_prior_beta_mix(weights, cooc , corpora, seq, alpha_prior, beta_prior, sampling = True):
    """This function is used to determine values of the prior mixture distribution.
    In the prior mixture distribution, individual components are Beta() distributions related to the probability 'p' of success: a mention of a compound (or specie) in an article also involving the MeSH descriptor 'M'  
    Weights used in the mixture model are probabilities that a mention between a compound and an article, that reached the targeted compound, involved the 'i' compound.

    Args:
        weights (list): A list of float used as weights in the mixture distribution.
        cooc (list): A list of integer values representing co-occurences between species in the graph and a particular MeSH descriptor  
        corpora (list):  A list of integer values representing copus sizes related to each compounds in the graph
        seq (float): the step used in np.arange to create the x vector of probabilities, only when sampling is True.
        alpha_prior (float): the alpha parameter of the beta prior distribution associated to the MeSH probabilities, p(M|S). Default to 1 for uninformative prior.
        beta_prior (float): the beta parameter of the beta prior distribution associated to the MeSH probabilities, p(M|S). Default to 1 for uninformative prior.
        sampling (bool, optional): Does the function have to compute a sampling of density values ?

    Returns:
        [collection]: A collection with:
            - alpha (list): vector of alpha parameters related to each individual Beta distribution in the prior mixture. 
            - beta (list): vector the beta parameters related to each individual Beta distribution in the prior mixture.
            - weights (list): vector of weights related to each individual Beta distribution in the prior mixture.
            - x (list): Probabilities
            - f (list): Densities 
            - mu (float): Mean of the distribution  
    """
    # Get parameters
    r = collections.namedtuple("priormix", ["alpha", "beta", "weights", "x", "f", "mu"])
    x = None
    y = None

    l = len(weights)

    # Get beta component paremeters for each compounds, using a posterior with uniformative prior
    alpha = [(cooc[it] + alpha_prior) for it in range(0, l)]
    beta = [(corpora[it] - cooc[it] + beta_prior) for it in range(0, l)]

    if sampling:
        x = np.arange(0, 1 + seq, seq).tolist()
        f_i = [ss.beta.pdf(x, a = alpha[it], b = beta[it]) for it in range(0, l)]
        # Create Beta mix:
        y = np.dot(weights, f_i)

    # Get mean
    mu_i = [(alpha[it]/(alpha[it] + beta[it])) for it in range(0, l)]
    mu = np.dot(weights, mu_i)

    mix = r(alpha, beta, weights, x, y, mu)

    return mix
    


def create_posterior_beta_mix(k, n, weights_pior, alpha_prior, beta_prior, seq, use_log = True, sampling = True):
    """This function is used to compute the posterior mixture distribution. The prior mixture model is updated for the observation of the coocurences on the targeted specie

    Args:
        k (integer): The coocurence for the MeSH 'M' for the targeted specie (number of success)
        n (integer): The corpus size for the targeted specie (number of trials)
        weights_pior (list): The prior mixture distribution weights (item 'weights' in create_prior_beta_mix result)
        alpha_prior (list): The prior mixture distribution alpha parameters (item 'alpha' in create_prior_beta_mix result)
        beta_prior (list): The prior mixture distribution beta parameters (item 'beta' in create_prior_beta_mix result)
        seq (float): The step used in np.arange to create the x vector of probabilities, only when sampling is True.
        use_log (boolean, optional): A boolean telling if the computation of the weights W have to be achieve using classic formula or using logs (As alpha and beta values are often large, log method is prefered). Default = True
        sampling (bool, optional): Does the function have to compute a sampling of density values ?

    Returns:
        [collection]: A collection with:
            - alpha (list): vector of alpha parameters related to each individual Beta distribution in the posterior mixture. 
            - beta (list): vector the beta parameters related to each individual Beta distribution in the posterior mixture.
            - weights (list): vector of weights related to each individual Beta distribution in the posterior mixture.
            - x (list): Probabilities
            - f (list): Densities 
            - mu (float): Mean of the distribution 
    """
    r = collections.namedtuple("posteriormix", ["alpha", "beta", "weights", "x", "f", "mu"])
    l = len(weights_pior)
    x = None
    y = None

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

    if sampling:
        x = np.arange(0, 1 + seq, seq).tolist()
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


def create_vizu_data(index, probabilities, q, weights, data, cptm = "c"):
    # Step 1: Get all potential contributors : repeat procedure for computing weights
    _sigmas = probabilities.to_numpy()
    _constrains = copy.copy(_sigmas)
    np.fill_diagonal(_constrains, 0)
    _constrains = _constrains @ np.diag(1/(_constrains.sum(axis=0)))
    potential_contributors = _constrains[index, :] > q
    # We also select the targeted compound:
    potential_contributors[index] = True
    # Get data
    w_data = data.loc[potential_contributors, ]
    # Force integer type before pasting data
    w_data = w_data.astype({"TOTAL_PMID_SPECIE": int, "COOC": int})
    # Add compartiment to specie for mapping: 
    w_data["SPECIE"] = w_data["SPECIE"] + "_" + cptm
    w_data["STATS"] = w_data["SPECIE"] + " (" + w_data["COOC"].astype(str) + "/" + w_data["TOTAL_PMID_SPECIE"].astype(str) + ")"
    w_data = w_data.drop(["index", "CID", "TOTAL_PMID_SPECIE", "COOC"], axis = 1)
    w_data = w_data.rename(columns={"SPECIE": "metaboliteDBIdentifier"})
    w_data.loc[index]
    return w_data

def compute_contributors_number(weights, labels):
    """This function is used to return the number of contributors

    Args:
        weights (np.array): the wegiht matrix
        labels (list): list of species labels

    Returns:
        (pd.Dataframe): a DataFrame containing for each specie the number of contributors
    """
    N = np.sum((weights > 0), axis = 0, dtype = np.int)
    res = pd.DataFrame({"SPECIE": labels, "NbCtb": N})
    
    return res

def compute_contributors_corpora_sizes(weights, table_species_corpora, labels):
    """
    This function is used to compute the weighted average corpora size of the contributors for each compounds, using weights from the weight matrix
    Args:
        weights (np.array): the wegiht matrix
        table_species_corpora (pd.Dataframe): The table containing corpora sizes
        labels (list): list of species labels
    
    Returns:
        (pd.Dataframe): a DataFrame containing for each specie the average corpora size of the contributors, with NaN if there is no available contributors for the compound
    """
    corpora = np.array([table_species_corpora["TOTAL_PMID_SPECIE"]]).T
    C = weights.T @ corpora
    # If all the weight are null, set NaN. Note that weights are necessarily null if corpora sizes for contributors are null
    C[np.sum(weights, axis = 0) == 0] = np.NaN
    C = np.round(C, decimals = 3)
    res = pd.DataFrame({"SPECIE": labels, "CtbAvgCorporaSize": C[:, 0]})
    
    return res

def compute_contributors_distances(weights, g, labels):
    """
    This function is used to compute the weighted average distance of the contributors for each compounds, using weights from the weight matrix
    Args:
        weights (np.array): the wegiht matrix
        g (igraph.Graph): the compound graph
        labels (list): list of species labels

    Returns:
        (pd.Dataframe): a DataFrame containing for each specie the average distance of the contributors, with NaN if there is no available contributors for the compound
    """
    D = g.shortest_paths()
    # By multiplying the weight matrix and the distance matrix element-wise, we can compute the values of the weighted average
    m = weights * D
    # To determine the average mean, we just need to sum these values
    M = m.sum(axis = 0)
    # The only way for the weighted average to be null is when all weights are null. In this case we set NaN
    M[M == 0] = np.NaN
    M = np.round(M, decimals = 3)
    res = pd.DataFrame({"SPECIE": labels, "CtbAvgDistance": M})

    return res


def E(w):
    """
    Compute Shanon Entropy
    Args:
        w ([list]): a vector of probabilities

    Returns:
        [float]: Entropy
    """
    if np.sum(w):
        # If there is only one contributor, the entropy will be computed as -log2(1) and will return -0.0 in numpy. To fix this, we return the abs value of the entropy
        return abs(-sum([ np.log2(p) * p for p in w if p != 0 ]))
    return np.NaN

def compute_Entropy_matrix(weight_matrix, labels):
    """
    This function is used to compute the Entropy values for each compound based on the contributor's distribution. If there is no contributor, the value returned is NaN.
    Args:
        weight_matrix (np.array): the weight matrix
        labels (list): list of species labels
    Returns:
        [list]: a list containing the entropy associated to the distribution of contributors for each compound
    """
    # When using list comprehension, python iter by row so we transpose the weight matrix to iter by columns
    Entropy = [E(w) for w in weight_matrix.T]
    Entropy = np.round(Entropy, decimals = 3)
    res = pd.DataFrame({"SPECIE": labels, "Entropy":Entropy})
    return res

def plot_distributions(prior_mix, posterior_mix):
    """This function is used to plot prior distribution against a posterior distribution

    Args:
        prior_mix (collection): A collection containing information about a prior distribution with x probabilties and associated densities from create_prior_beta_mix or simple_prior with the samping = True
        posterior_mix (collection): A collection containing information about a posterior distribution with x probabilties and associated densities from create_posterior_beta_mix or simple_posterior with the samping = True
    """
    plt.plot(prior_mix.x, prior_mix.f, label = "prior", color = "blue")
    plt.plot(posterior_mix.x, posterior_mix.f, label = "Posterior", color = "red")
    plt.title('Differences between prior mix and posterior mix distribution ')
    plt.legend()
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.show()

def plot_prior_mix_distributions(prior_mix, labels, seq, top = 10):
    """This function is used to plot distribution of each components of a prior mixture distribution. The function compute itself the densities of each component. Since there could be dozens of contributors, we only plot the top n (top argument) for clarity.

    Args:
        prior_mix (collections): A collection containing information about the  prior mixture distribution: weights, alpha and beta parameters
        labels (list): A list of compound (or specie) labels associated to each component of the mixture. 
        seq (float): The step used in np.arange to create the x vector of probabilities.
        top (int): The top n (maximum) of contributors that should be plotted 
    """
    x = np.arange(0, 1 + seq, seq).tolist()
    weights = prior_mix.weights
    # To plot only the top n contributors, we first order the index of weights in decreasing order
    ordered_i_w = np.flip(np.argsort(weights))
    # We go trough the list of weight until we reach the n'th contributor, or the last contributor if there are les than n
    for i in range(0, min(len(weights), top)):
        it = ordered_i_w[i]
        f = ss.beta.pdf(x, a = prior_mix.alpha[it], b = prior_mix.beta[it])
        y = weights[it] * f
        plt.plot(x, y, label = labels[it] + ": Beta(" + str(round(prior_mix.alpha[it], 2)) + ", " + str(round(prior_mix.beta[it], 2)) + ") - w = " + str(round(weights[it],2)))
    plt.title("Top " + str(min(len(weights), top)) + " - Prior decomposition")
    plt.figtext(0.995, 0.01, 'w is the weight of the component represented by the species in the beta mixture distribution', ha='right', va='bottom')
    plt.legend()
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.show()


def compute_mix_CDF(p, weights, alpha, beta):
    """This function is used to compute the CDF of a mixture distribution

    Args:
        p (float): the probability P(x <= p)
        weights (list): The mixture distribution weights
        alpha (list): The mixture distribution alpha parameters
        beta (list): The mixture distribution beta parameters

    Returns:
        [type]: [description]
    """
    cdf_i = [ss.beta.cdf(p, alpha[it], beta[it]) for it in range(0, len(weights))]
    cdf = np.dot(weights, cdf_i)
    return cdf

def computation(index, data, p, alpha_prior, beta_prior, seq = 0.0001, plot = False, weigth_limit = 1e-5):
    """This function is used to compute the complete analysis for a Compound - MeSH relation.
    If the neighborhood can't provide information about the prior distribution, then the default prior from estimate_prior_distribution_mesh_V2 is used, otherwise we will used the prior mixture.

    Args:
        index (integer): the index of the specie if the metabolic network
        data (pandas.DataFrame): data related to corpus size of each compound in the metabolic network and their co-occurence with the studied MeSH 
        p (float): The general probability to observed a mention of a compound in an article, also involving the MeSH.
        alpha_prior (float): The alpha parameter of the MeSH's prior distribution (Cf. estimate_prior_distribution_mesh_V2)
        beta_prior (float): The beta parameter of the MeSH's prior distribution (Cf. estimate_prior_distribution_mesh_V2)
        seq (float, optional): The step used to create a x vector of probabilities (used for plotting distribution only). Defaults to 0.0001.
        plot (bool, optional): Does the function has to plot prior and posterior distributions ?. See plot_prior_mix_distributions and plot_distributions. Defaults to False.
        weigth_limit (float, optional): If the weight of a compound in the prior mixture is lower than this threshild, the compound is removed from the mixture. It may be usefull when plotting distribution as there could be a lot of compounds involved in the mxiture. Defaults to 1e-5.

    Returns:
        [collection]: A collection with:
        - Mean (float): The mean of the posterior distribution. 
        - CDF (float): The probability P(q <= p(M)) derived from the CDF of the posterior distribution. The more this probability is low, the more we are certain that the mean of the posterior distribution is higher than the general probability to observed the MeSH (the 'p' argument of the function), representing independence hypothsis.
        - Log2FC (float): The log2 fold change between the mean of the posterior distribution and the general probability to observed the MeSH (the 'p' argument of the function)
        - priorCDF: Same as CDF, but for the mixture prior.
        - priorLog2FC: Same as Log2FC, but for the mixture prior.
    """
    
    # Out
    r = collections.namedtuple("out", ["Mean", "CDF", "Log2FC", "priorCDF", "priorLog2FC", "NeighborhoodInformation"])

    weights = data["weights"].tolist()
    del weights[index]
    cooc = data["COOC"].tolist()
    k = cooc.pop(index)
    corpora = data["TOTAL_PMID_SPECIE"].tolist()
    labels = data["SPECIE"].to_list()
    del labels[index]
    n = corpora.pop(index)
    # If all weights are null, no neighborhood information:
    if sum(weights) == 0:
        prior = simple_prior(alpha_prior, beta_prior, seq, sampling = plot)
        posterior = simple_posterior(k, n, alpha_prior, beta_prior, seq, sampling = plot)
        # In case of no neighborhood information, we simply plot prior vs posterior distributions:
        if plot:
            plot_distributions(prior, posterior)
        # Compute additional values:
        Log2numFC = np.log2(posterior.mu/p)
        cdf_posterior = ss.beta.cdf(p, posterior.alpha, posterior.beta)
        resultat = r(posterior.mu, cdf_posterior, Log2numFC, np.NaN, np.NaN, False)
        return resultat

    # Null weights have to be removed before the computation as we cill the log(weights) during the computation.
    to_remove = list()
    for it in range(0, len(weights)):
        # Check weight
        if not weights[it]:
            to_remove.append(it)
    # Once the list of items to be deleted is complete, we remove them. We need to delete them in reverse order so that we don't throw off the subsequent indexes.
    for rmv in sorted(to_remove, reverse=True):
        del weights[rmv]
        del cooc[rmv]
        del corpora[rmv]
        del labels[rmv]
    
    # Use initial prior on MeSH (uninformative or from glm) to build a prior mix using neighboors' observations
    prior_mix = create_prior_beta_mix(weights, cooc, corpora, seq, alpha_prior, beta_prior, sampling = plot)
    # Get ratio between initial prior on MeSH and (posterior) prior using neighboors' indicating whether the neighbours are in favour of the relationship
    prior_mix_CDF = compute_mix_CDF(p, prior_mix.weights, prior_mix.alpha, prior_mix.beta)
    
    # Compute prior Mean ratio
    prior_mean_ratio = np.log2(prior_mix.mu/p)
    # Posterior mix:
    posterior_mix = create_posterior_beta_mix(k, n, prior_mix.weights, prior_mix.alpha, prior_mix.beta, seq, sampling = plot)
    cdf_posterior_mix = compute_mix_CDF(p, posterior_mix.weights, posterior_mix.alpha, posterior_mix.beta)
    # print(labels)
    # print(prior_mix.weights)
    # print(prior_mix.alpha)
    # print(prior_mix.beta)
    # print("========================")
    # print(posterior_mix.weights)
    # print(posterior_mix.alpha)
    # print(posterior_mix.beta)

    Log2numFC = np.log2(posterior_mix.mu/p)

    if plot: 
        plot_prior_mix_distributions(prior_mix, labels, seq)
        plot_distributions(prior_mix, posterior_mix)
    
    resultat = r(posterior_mix.mu, cdf_posterior_mix, Log2numFC, prior_mix_CDF, prior_mean_ratio, True)

    return resultat


def specie_mesh(index, table_cooc, table_species_corpora, weights, table_mesh, forget):
    """This function is used to computed associations from a specific specie against all available MeSHs.

    Args:
        index (int): index of the specie in the metabolic network
        table_cooc (pandas.DataFrame): table of co-occurences
        table_species_corpora (pandas.DataFrame): table of specie corpora
        weights (numpy): weight matrix
        table_mesh (pandas.DataFrame): table of MeSH corpora
        forget (bool): Keep only prior information from the neighborhood, removing specie's observation
    """
    # Create result Dataframe from MeSH list
    mesh_list = table_mesh["MESH"].tolist()
    indexes = range(0, len(mesh_list))
    df_ = pd.DataFrame(index = indexes, columns = ["Mean", "CDF", "Log2FC", "priorCDF", "priorLog2FC", "NeighborhoodInformation"])

    # Prepare data table
    table_species_corpora["weights"] = weights[:, index].tolist()
    
    with progressbar.ProgressBar(max_value=len(indexes)) as bar:
        for i in indexes:
            mesh = mesh_list[i]
            # Get cooc vector. It only contains species that have at least one article, need to left join.
            cooc = table_cooc[table_cooc["MESH"] == mesh][["index", "COOC"]]
            # Get data
            data = pd.merge(table_species_corpora, cooc, on = "index", how = "left").fillna(0)
            
            # If forget option is true, remove observation of the studied specie
            if forget:
                data.loc[data["index"] == index, ["TOTAL_PMID_SPECIE", "COOC"]] = [0, 0]
            
            # Get MeSH info
            MeSH_info = table_mesh[table_mesh["MESH"] == mesh]
            p = float(MeSH_info["P"])
            
            # Computation
            r = computation(index, data, p, float(MeSH_info["alpha_prior"]), float(MeSH_info["beta_prior"]), seq = 0.0001)
            df_.loc[i] = r
            bar.update(i)
    
    df_.insert(0, "MESH", mesh_list)
    return(df_)

def mesh_specie(mesh, table_cooc, table_species_corpora, weights, table_mesh, forget):
    """This function is used to computed associations from a specific MeSH against all available species.

    Args:
        index (int): index of the specie in the metabolic network
        table_cooc (pandas.DataFrame): table of co-occurences
        table_species_corpora (pandas.DataFrame): table of specie corpora
        weights (numpy): weight matrix
        table_mesh (pandas.DataFrame): table of MeSH corpora
        forget (bool): Keep only prior information from the neighborhood, removing specie's observation.
    """
    specie_list = table_species_corpora["SPECIE"].tolist()
    indexes = range(0, len(specie_list))
    df_ = pd.DataFrame(index = indexes, columns = ["Mean", "CDF", "Log2FC", "priorCDF", "priorLog2FC", "NeighborhoodInformation"])
    
    # Get MeSH info
    MeSH_info = table_mesh[table_mesh["MESH"] == mesh]
    cooc = table_cooc[table_cooc["MESH"] == mesh][["index", "COOC"]]
    p = float(MeSH_info["P"])
    
    # Browser all species
    with progressbar.ProgressBar(max_value=len(indexes)) as bar:
        for i in indexes:
            table_species_corpora["weights"] = weights[:, i].tolist()
            data = pd.merge(table_species_corpora, cooc, on = "index", how = "left").fillna(0)
            
            # If forget option is true, remove observation of the studied specie
            if forget:
                data.loc[data["index"] == i, ["TOTAL_PMID_SPECIE", "COOC"]] = [0, 0]
            
            # Computation
            r = computation(i, data, p, float(MeSH_info["alpha_prior"]), float(MeSH_info["beta_prior"]), seq = 0.0001)
            df_.loc[i] = r
            bar.update(i)
    
    df_.insert(0, "SPECIE", specie_list)
    return(df_)

def association_file(f, table_cooc, table_species_corpora, weights, table_mesh, forget):
    """This function is used to compute all associations specified in a Dataframe (SPECIE, MESH)

    Args:
        f (pandas.Dataframe): A two columns Dataframe storing all SPECIE - MESH pairs that need to be computed.
        table_cooc (pandas.DataFrame): table of co-occurences
        table_species_corpora (pandas.DataFrame): table of specie corpora
        weights (numpy): weight matrix
        table_mesh (pandas.DataFrame): table of MeSH corpora
        forget (bool): Keep only prior information from the neighborhood, removing specie's observation.

    Returns:
        [type]: [description]
    """
    associations = pd.concat([f, pd.DataFrame(columns = ["Mean", "CDF", "Log2FC", "priorCDF", "priorLog2FC", "NeighborhoodInformation"])])
    n = len(associations)
    
    # Browse associations
    with progressbar.ProgressBar(max_value = n) as bar:
        for i in range(0, n):
            specie = str(associations.iloc[[i], 0].item())
            mesh = str(associations.iloc[[i], 1].item())
            index = int(table_species_corpora[table_species_corpora["SPECIE"] == specie]["index"])
            # Prepare data
            table_species_corpora["weights"] = weights[:, index].tolist()
            cooc = table_cooc[table_cooc["MESH"] == mesh][["index", "COOC"]]
            data = pd.merge(table_species_corpora, cooc, on = "index", how = "left").fillna(0)
            
            # If forget option is true, remove observation of the studied specie
            if forget:
                data.loc[data["index"] == index, ["TOTAL_PMID_SPECIE", "COOC"]] = [0, 0]
            
            # Get MeSH info
            MeSH_info = table_mesh[table_mesh["MESH"] == mesh]
            p = float(MeSH_info["P"])

            # Computation
            r = computation(index, data, p, float(MeSH_info["alpha_prior"]), float(MeSH_info["beta_prior"]), seq = 0.0001, plot = False)
            associations.iloc[i, 2:8] = list(r)
            bar.update(i)
    return associations