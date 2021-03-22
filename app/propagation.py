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
    # Create a table to map species labels (SPECIE column) to indexes in the graph
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
        alpha (float): The damping factor.
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
    # print(str(c) + " iterations to convergence.")

    # Insert 0 at targeted index
    r = np.insert(new_pi, i, 0, axis = 0)
    # Float are basically imprecise and so after several matrix multiplications, the sum of probabilities in the vector may not equal to 1, but 0.9999999999999999 or 1.0000000000000001 for example. 
    if np.sum(r, axis = 0, dtype = np.float16) != 1:
        print("Warning at index " + str(i) + ": the final probability vector does not sum to 1. This may be due to float approximation errors")

    return r

def propagation_volume(g, alpha = 0.8, name_att = "label", direction = "both"):
    """This function is used to compute the PPR, excluding the targeted node itself, for each node of the graph

    Args:
        g (igraph.Graph): The compound graph
        alpha (float, optional): The damping factor. Defaults to 0.8.
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
        full[:, i] = compute_PR(A, i, alpha)[:, 0]
    
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
    """This function is used to estimate parameters and distribution of a simple prior distribution, no mixture.

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
    """This function is used to estimate parameters and distribution of a simple posterior distribution, no mixture.

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

def plot_prior_mix_distributions(prior_mix, labels, seq):
    """This function is used to plot distribution of each components of a prior mixture distribution. The function compute itself the densities of each component.

    Args:
        prior_mix (collections): A collection containing information about the  prior mixture distribution: weights, alpha and beta parameters
        labels (list): A list of compound (or specie) labels associated to each component of the mixture. 
        seq (float): The step used in np.arange to create the x vector of probabilities.
    """
    x = np.arange(0, 1 + seq, seq).tolist()
    weights = prior_mix.weights
    for it in range(0, len(weights)):
        f = ss.beta.pdf(x, a = prior_mix.alpha[it], b = prior_mix.beta[it])
        y = weights[it] * f
        plt.plot(x, y, label = labels[it])
    plt.title('Prior decomposition')
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
        alpha_prior (float): The alpha parameter of the prior distribution
        beta_prior (float): The beta parameter of the prior distribution
        seq (float, optional): The step used to create a x vector of probabilities (for plotting distribution only). Defaults to 0.0001.
        plot (bool, optional): Does the function has to plot prior and posterior distributions ?. See plot_prior_mix_distributions and plot_distributions. Defaults to False.
        weigth_limit (float, optional): If the weight of a compound in the prior mixture is lower than this threshild, the compound is removed from the mixture. it may be usefull when plotting distribution as there could be a lot of compounds involved in the mxiture. Defaults to 1e-5.

    Returns:
        [collection]: A collection with:
        - Mean (float): The mean of the posterior distribution. 
        - CDF (float): The probability P(p <= p(M)) derived from the CDF of the posterior distribution. The more this probability is low, the more we are certain that the mean of the posterior distribution is higher than the general probability to observed the MeSH (independance hypothesis)
        - Log2FC (float): The log2 fold change between the mean of the posterior distribution and the general probability to observed the MeSH
        - priorCDFratio: The log2 ratio of the CDF probabilities P(p <= p(M)) obtained between the initial prior and the mixture prior. When this value is high, it indicates that the studied MeSH is more frequent than usual in the neiborhood of the targeted compound. This value is correlated with the CDF. This value is NaN is the neighborhood can't provide information, as these both prior will be the same
    """
    
    # Out
    r = collections.namedtuple("out", ["Mean", "CDF", "Log2FC", "priorCDFratio"])

    weights = data["weights"].tolist()
    del weights[index]
    cooc = data["COOC"].tolist()
    k = cooc.pop(index)
    corpora = data["TOTAL_PMID_SPECIE"].tolist()
    labels = data["SPECIE"].to_list()
    del labels[index]
    n = corpora.pop(index)

    # If all weights are null, no neighborhood information:
    if all(w == 0 for w in weights):
        print("Neiborhood literature information does not reach the targeted compound. You should increase the damping factor. Use default prior.")
        prior = simple_prior(alpha_prior, beta_prior, seq, sampling = plot)
        posterior = simple_posterior(k, n, alpha_prior, beta_prior, seq, sampling = plot)
        # In case of no neighborhood information, we simply plot prior vs posterior distributions:
        if plot:
            plot_distributions(prior, posterior)
        # Compute additional values:
        Log2numFC = np.log2(posterior.mu/p)

        resultat = r(posterior.mu, ss.beta.cdf(p, posterior.alpha, posterior.beta), Log2numFC, np.NaN)
        return resultat

    # Check for other null weights, in case of low alpha (damping) for instance. We consider a weight null if weight < 1e-5
    if not all(w > weigth_limit for w in weights):
        to_remove = list()
        for it in range(0, len(weights)):
            # Check if weights == 0
            if weights[it] <= weigth_limit:
                # print("Warning: weight at index " + str(it) + " is null: " + labels[it] + " Low damping factor used ?")
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
    # If prior mix CDF is already estimated to 0, set log2FC to infinite
    if not prior_mix_CDF:
        print("Warning: prior mix CDF is estimated to 0. The value of the CDF ratio between MeSH prior and prior mix is set to Inf.")
        prior_cdf_ratios = np.Inf
    else:
        prior_cdf_ratios = np.log2(ss.beta.cdf(p, alpha_prior, beta_prior)/prior_mix_CDF)
    # Posterior mix:
    posterior_mix = create_posterior_beta_mix(k, n, prior_mix.weights, prior_mix.alpha, prior_mix.beta, seq, sampling = plot)
    cdf_posterior_mix = compute_mix_CDF(p, posterior_mix.weights, posterior_mix.alpha, posterior_mix.beta)
    # print(labels)
    # print(prior_mix.weights)
    # print(prior_mix.alpha)
    # print(prior_mix.beta)
    # print("========================")
    # print(posterior_mix.weights)
    # print(posterior_mix.alpha)
    # print(posterior_mix.beta)

    Log2numFC = np.log2(posterior_mix.mu/p)
    
    if plot: 
        plot_prior_mix_distributions(prior_mix, labels, seq)
        plot_distributions(prior_mix, posterior_mix)
    
    resultat = r(posterior_mix.mu, cdf_posterior_mix, Log2numFC, prior_cdf_ratios)

    return resultat


def specie_mesh(index, table_cooc, table_corpora, weights, table_mesh):
    
    print("Treat index: " + str(index))
    # Create result Dataframe from MeSH list
    mesh_list = table_mesh["MESH"].tolist()
    indexes = range(0, len(mesh_list))
    df_ = pd.DataFrame(index = indexes, columns = ["Mean", "CDF", "Log2FC", "priorCDFratio"])

    # Prepare data table
    table_corpora.insert(2, "weights", weights[:, index].tolist())
    
    for i in indexes:
        print(str(i) + "/" + str(len(indexes)))
        mesh = mesh_list[i]
        # Get cooc vector. It only contains species that have at least one article, need to left join.
        cooc = table_cooc[table_cooc["MESH"] == mesh][["index", "COOC"]]
        # Get data
        data = pd.merge(table_corpora, cooc, on = "index", how = "left").fillna(0)
        MeSH_info = table_mesh[table_mesh["MESH"] == mesh]
        p = float(MeSH_info["P"])
        r = computation(index, data, p, float(MeSH_info["alpha_prior"]), float(MeSH_info["beta_prior"]), seq = 0.0001)
        df_.loc[i] = r
    
    df_.insert(0, "MESH", mesh_list)
    return(df_)
        

# if data["COOC"][index] > 0:
# Test of interest (mean raw Log2FC neighboors > 1):
# testFC = np.dot(data[data["TOTAL_PMID_SPECIE"] != 0]["weights"], ((data[data["TOTAL_PMID_SPECIE"] != 0]["COOC"]/data[data["TOTAL_PMID_SPECIE"] != 0]["TOTAL_PMID_SPECIE"])/p))
# if testFC != 0 and np.log2(testFC) > 1 :
# if r.CDF <= 0.001:



