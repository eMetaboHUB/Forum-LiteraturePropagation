import sys, os
import igraph as ig
import pandas as pd
import numpy as np
import collections
import warnings
import scipy.special as sc
import scipy.stats as ss
import progressbar
from matplotlib import cm
np.set_printoptions(suppress=True)

from plots import *


def observation_uninformative_prior(k, n, seq, sampling=True):
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
        y = ss.beta.pdf(x, a=alpha, b=beta)

    # Get mean of the distribution (Bayes estimator of p)
    mu = (k + 1)/(n + 2)

    res = r(alpha, beta, x, y, mu)

    return res

def estimate_prior_distribution_mesh(mu, sample_size):
    """This function is used to estimate the prior distribution of the probability that a mention of a compound discuss the studied MeSH.
    This prior distribution is based on the independance hypothesis.
    Args:
        mu (float): The expected mean of the prior probability distributions
        sample_size (integer): The expected sample size from which this prior distribution has been generated, related to the variability, certainity, of the distribution. This parameter could also be refered as the concentration parameter. 

    Returns:
        [collection]: A collection with:
            - alpha: the alpha parameter of the prior distribution
            - beta: the beta parameter of the prior distribution
    """
    r = collections.namedtuple("prior_mesh", ["alpha", "beta"])

    alpha = mu * sample_size
    beta = (1 - mu) * sample_size

    result = r(alpha, beta)
    return result

def simple_prior(alpha_prior, beta_prior, seq, sampling=True):
    """This function is used to compute the density of a simple prior distribution, no mixture.

    Args:
        alpha_prior (float): The alpha parameter of the prior probability distribution (Cf. estimate_prior_distribution_mesh)
        beta_prior (float): The beta parameter of the prior probability distribution (Cf. estimate_prior_distribution_mesh)
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
        y = ss.beta.pdf(x, a=alpha_prior, b=beta_prior)

    # Get mean
    mu = alpha_prior/(alpha_prior + beta_prior)
    res = r(alpha_prior, beta_prior, x, y, mu)

    return res

def simple_posterior(cooc, corpora, alpha_prior, beta_prior, seq, sampling=True):
    """This function is used to estimate parameters and density of a simple posterior distribution, no mixture.

    Args:
        cooc (list): The co-occurence between the specie in the graph and a particular MeSH descriptor
        corpora (list): The corpus size related to the specie in the graph
        alpha_prior (float): The alpha parameter of the prior probability distribution (Cf. estimate_prior_distribution_mesh)
        beta_prior (float): The beta parameter of the prior probability distribution (Cf. estimate_prior_distribution_mesh)
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
        y = ss.beta.pdf(x, a=alpha, b=beta)

    # Get mean
    mu = alpha/(alpha + beta)
    res = r(alpha, beta, x, y, mu)

    return res

def create_prior_beta_mix(weights, cooc, corpora, seq, alpha_prior, beta_prior, sampling=True):
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
    r = collections.namedtuple("priormix", ["alpha", "beta", "weights", "x", "f", "mu", "l_mu"])
    x = None
    y = None

    l = len(weights)

    # Get beta component paremeters for each compounds, using a posterior with uniformative prior
    alpha = [(cooc[it] + alpha_prior) for it in range(0, l)]
    beta = [(corpora[it] - cooc[it] + beta_prior) for it in range(0, l)]

    if sampling:
        x = np.arange(0, 1 + seq, seq).tolist()
        f_i = [ss.beta.pdf(x, a=alpha[it], b=beta[it]) for it in range(0, l)]

        # Create Beta mix:
        y = np.dot(weights, f_i)

    # Get mean
    mu_i = [(alpha[it]/(alpha[it] + beta[it])) for it in range(0, l)]
    mu = np.dot(weights, mu_i)

    mix = r(alpha, beta, weights, x, y, mu, mu_i)

    return mix


def create_posterior_beta_mix(k, n, weights_pior, alpha_prior, beta_prior, seq, use_log=True, sampling=True):
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
    r = collections.namedtuple("posteriormix", ["alpha", "beta", "weights", "x", "f", "mu", "l_mu"])
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
        W = [round((np.exp((Z[it] - (sc.logsumexp(Z))))), 9) for it in range(0, l)]

    else:
        C = [sc.beta(alpha_post[it], beta_post[it])/sc.beta(alpha_prior[it], beta_prior[it]) for it in range(0, l)]
        W = [round((weights_pior[it] * C[it]/(np.dot(weights_pior, C))), 9) for it in range(0, l)]

    if sampling:
        x = np.arange(0, 1 + seq, seq).tolist()

        # Create posterior distribution by componennts
        f_post_i = [ss.beta.pdf(x, a=alpha_post[it], b=beta_post[it]) for it in range(0, l)]

        # Get posteriors probabilities
        y = np.dot(W, f_post_i)

    # Get mean
    mu_i = [(alpha_post[it]/(alpha_post[it] + beta_post[it])) for it in range(0, l)]
    mu = np.dot(W, mu_i)

    mix = r(alpha_post, beta_post, W, x, y, mu, mu_i)

    return mix

def compute_contributors_number(weights, labels):
    """This function is used to return the number of contributors

    Args:
        weights (np.array): the wegiht matrix
        labels (list): list of species labels

    Returns:
        (pd.Dataframe): a DataFrame containing for each specie the number of contributors
    """

    N = np.sum((weights > 0), axis=0, dtype=np.int)
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
    C[np.sum(weights, axis=0) == 0] = np.NaN
    C = np.round(C, decimals=3)
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
    M = m.sum(axis=0)

    # The only way for the weighted average to be null is when all weights are null. In this case we set NaN
    M[M == 0] = np.NaN
    M = np.round(M, decimals=3)
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
        return abs(-sum([np.log2(p) * p for p in w if p != 0]))

    return np.NaN

def compute_entropy_matrix(weight_matrix, labels):
    """
    This function is used to compute the Entropy values for each compound based on the contributor's distribution. If there is no contributor, the value returned is NaN.
    Args:
        weight_matrix (np.array): the weight matrix
        labels (list): list of species labels
    Returns:
        [list]: a list containing the entropy associated to the distribution of contributors for each compound
    """

    # When using list comprehension, python iter by row so we transpose the weight matrix to iter by columns
    entropy = [E(w) for w in weight_matrix.T]
    entropy = np.round(entropy, decimals=3)
    res = pd.DataFrame({"SPECIE": labels, "Entropy":entropy})

    return res

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

def compute_log_odds(cdf):
    """
    This function is used to compute the log(odds) from the CDF.
    The success probability is defined as (1-CDF) so P(p > mu)
    If CDF = 0, reports infinite odds.
    If CDF = 1, reports -infinite odds
    Args:
        cdf (float): The computed CDF from *computation* (P(p <= mu))

    Returns:
        [type]: [description]
    """
    # Due to float approximation, while computing the CDF of the mixture, it could appears > 1 due to approximation. To avoid issues with log next, we set it to 1. Rounding would not help.
    if cdf > 1:
        cdf = 1

    with np.errstate(all='ignore'):
        log_odds = np.log((1 - cdf)) - np.log(cdf)

    return log_odds



def computation(index, data, p, alpha_prior, beta_prior, seq=0.0001, report=None, weigth_limit=1e-5, species_name_path=None, update_data=False):
    """This function is used to compute the complete analysis for a Compound - MeSH relation.
    If the neighborhood can't provide information about the prior distribution, then the default prior from estimate_prior_distribution_mesh is used, otherwise we will used the prior mixture.

    Args:
        index (integer): the index of the specie if the metabolic network
        data (pandas.DataFrame): data related to corpus size of each compound in the metabolic network and their co-occurence with the studied MeSH
        p (float): The general probability to observed a mention of a compound in an article, also involving the MeSH.
        alpha_prior (float): The alpha parameter of the MeSH's prior distribution (Cf. estimate_prior_distribution_mesh)
        beta_prior (float): The beta parameter of the MeSH's prior distribution (Cf. estimate_prior_distribution_mesh)
        seq (float, optional): The step used to create a x vector of probabilities (used for plotting distribution only). Defaults to 0.0001.
        report (optional): Path to output the html report. Defaults to None.
        weigth_limit (float, optional): If the weight of a compound in the prior mixture is lower than this threshild, the compound is removed from the mixture. It may be usefull when plotting distribution as there could be a lot of compounds involved in the mxiture. Defaults to 1e-5.
        species_name_path (str): Path to the file containing species' names for figures legend
        update_data (bool, optional): Does the data table need to be updated with posterior weights, cdf, etc of each contributors (for export)

    Returns:
        [dict]: A dictionnary with:
        - TOTALspecies_name_path_PMID_SPECIE (int): The total number of mentions for the targeted compound
        - COOC (int): The total number of co-occurences between the targeted compound and the MeSH
        - Mean (float): The mean of the posterior distribution.
        - CDF (float): The probability P(q <= p(M)) derived from the CDF of the posterior distribution. The more this probability is low, the more we are certain that the mean of the posterior distribution is higher than the general probability to observed the MeSH (the 'p' argument of the function), representing independence hypothsis.
        - LogOdds (float): The logarithm of the Odds, computed from the CDF
        - Log2FC (float): The log2 fold change between the mean of the posterior distribution and the general probability to observed the MeSH (the 'p' argument of the function)
        - priorLogOdds: Same as LogOdds, but for the mixture prior.
        - priorLog2FC: Same as Log2FC, but for the mixture prior.
    """

    # Get data and remove the line corresponding to the targeted specie in data
    k = int(data.loc[index, "COOC"])
    n = int(data.loc[index, "TOTAL_PMID_SPECIE"])
    data.drop(index=index, inplace=True)

    # If all weights are null, no neighborhood information.
    if sum(data["weights"]) == 0:

        # use initial prior
        prior = simple_prior(alpha_prior, beta_prior, seq, sampling=report)
        log2fc_prior = np.log2(prior.mu/p)
        cdf_prior = ss.beta.cdf(p, prior.alpha, prior.beta)
        log_odds_prior = compute_log_odds(cdf_prior)

        # If there are observations available:
        if n > 0:
            posterior = simple_posterior(k, n, alpha_prior, beta_prior, seq, sampling = report)
            # Compute Log2FC:
            log2fc = np.log2(posterior.mu/p)

            # Compute the CDF from the posterior distribution
            cdf_posterior = ss.beta.cdf(p, posterior.alpha, posterior.beta)

            # Compute Log(odds) from the CDF
            log_odds = compute_log_odds(cdf_posterior)

            resultat = dict(zip(["TOTAL_PMID_SPECIE", "COOC", "Mean", "CDF", "LogOdds", "Log2FC", "priorLogOdds", "priorLog2FC", "NeighborhoodInformation"], [n, k, posterior.mu, cdf_posterior, log_odds, log2fc, log_odds_prior, log2fc_prior, False]))

            # In case of no neighborhood information, we simply plot prior vs posterior distributions:
            if report:
                f1 = plot_distributions_plotly(prior, posterior)
                generate_html_report(report, [f1], ["Prior .VS. Posterior"], resultat, pd.DataFrame())

        # If there are no observations:
        else:

            resultat = dict(zip(["TOTAL_PMID_SPECIE", "COOC", "Mean", "CDF", "LogOdds", "Log2FC", "priorLogOdds", "priorLog2FC", "NeighborhoodInformation"], [n, k, prior.mu, cdf_prior, log_odds_prior, log2fc_prior, np.NaN, np.NaN, False]))

            if report:
                f1 = go.Figure()
                f1.add_trace(go.Scatter(x=prior.x, y=prior.f, line=dict(color="blue"), name="Prior"))
                f1.update_layout(title="Prior from the MeSH overall frequency", xaxis=dict(title="Probability", titlefont_size=25, tickfont_size=20), yaxis=dict(title="Density", titlefont_size=25, tickfont_size=20), template="simple_white")
                generate_html_report(report, [f1], ["MeSH overall prior"], resultat,  pd.DataFrame())

        return resultat

    # Null weights have to be removed before the computation as we will use the log(weights) during the computation.
    to_remove = data[data["weights"] == 0].index
    data.drop(index=to_remove, inplace=True)
    data.reset_index(drop=True, inplace=True)

    # Use initial prior on MeSH (uninformative or from glm) to build a prior mix using neighboors' observations
    prior_mix = create_prior_beta_mix(data["weights"].tolist(), data["COOC"].tolist(), data["TOTAL_PMID_SPECIE"].tolist(), seq, alpha_prior, beta_prior, sampling=report)

    # Get ratio between initial prior on MeSH and (posterior) prior using neighboors' indicating whether the neighbours are in favour of the relationship
    prior_mix_CDF = compute_mix_CDF(p, prior_mix.weights, prior_mix.alpha, prior_mix.beta)

    # Compute priorLogOdds
    prior_log_odds = compute_log_odds(prior_mix_CDF)

    # Compute Log2FC ratio from prior mixture
    prior_log2fc = np.log2(prior_mix.mu/p)

    # If there are observations available:
    if n > 0:

        # Posterior mix:
        posterior_mix = create_posterior_beta_mix(k, n, prior_mix.weights, prior_mix.alpha, prior_mix.beta, seq, sampling=report)

        # Compute CDF of the posterior distibution
        cdf_posterior_mix = compute_mix_CDF(p, posterior_mix.weights, posterior_mix.alpha, posterior_mix.beta)

        # Compute Log(odds) from the CDF
        log_odds = compute_log_odds(cdf_posterior_mix)

        # Compute log2fc from the posterior mixture
        log2fc = np.log2(posterior_mix.mu/p)

        # If the inputs are a specific specie with a specific MeSH, we need to export the data table
        if update_data:

            # Rename weight column
            data.rename(columns={'weights':'PriorWeights'}, inplace=True)

            # As null weight have been removed during computation, we use SPECIE instead of index as key
            data["PostWeights"] = float(0)
            data["PriorLogOdds"] = np.NaN
            data["PostLogOdds"] = np.NaN
            data["PriorLog2FC"] = np.NaN
            data["PostLog2FC"] = np.NaN
            for j in data.index:
                # Contributor posterior weight
                data.loc[j, "PostWeights"] = posterior_mix.weights[j]

                # Contributor prior LogLodds
                ctb_prior_cdf = ss.beta.cdf(p, prior_mix.alpha[j], prior_mix.beta[j])
                data.loc[j, "PriorLogOdds"] = compute_log_odds(ctb_prior_cdf)

                # Contributor posterior LogLodds
                ctb_post_cdf = ss.beta.cdf(p, posterior_mix.alpha[j], posterior_mix.beta[j])
                data.loc[j, "PostLogOdds"] = compute_log_odds(ctb_post_cdf)

                # Contributor prior Log2FC
                data.loc[j, "PriorLog2FC"] = np.log2(prior_mix.l_mu[j]/p)

                # Contributor posterior LogFC
                data.loc[j, "PostLog2FC"] = np.log2(posterior_mix.l_mu[j]/p)

        resultat = dict(zip(["TOTAL_PMID_SPECIE", "COOC", "Mean", "CDF", "LogOdds", "Log2FC", "priorLogOdds", "priorLog2FC", "NeighborhoodInformation"], [n, k, posterior_mix.mu, cdf_posterior_mix, log_odds, log2fc, prior_log_odds, prior_log2fc, True]))

        # Plot figure ? Only when the inputs are a specific specie with a specific MeSH
        if report:

            # If names have been provided, use them instead of species labels in Figures:
            if "SPECIE_NAME" in data.columns:
                names = data["SPECIE_NAME"].tolist()
            else:
                names = data["SPECIE"]


            # We set the top to top 10:
            top = 10

            # We select the union of the top 10 contributors in the both groups and then assign a unique color to it in a dict (so max number of contributors displayed is 20)
            set_contributors = set([names[i] for i in np.argsort(prior_mix.weights)[::-1][:top]] + [names[i] for i in np.argsort(posterior_mix.weights)[::-1][:top]])

            # We need to keep the same color palette between the both plots
            palette = dict(zip(set_contributors, cm.tab20(np.linspace(0, 1, len(set_contributors)))))

            # Plot
            # f1 = plot_mix_distributions_plotly(prior_mix, names, seq, "Prior components", palette, top)
            # f2 = plot_mix_distributions_plotly(posterior_mix, names, seq, "Posterior components", palette, top)
            # f3 = plot_distributions_plotly(prior_mix, posterior_mix)

            # Contribution plots
            f4 = contributions_plot(data, names, "PriorWeights", "PriorLogOdds")
            f5 = contributions_plot(data, names, "PostWeights", "PostLogOdds")

            # Generate report:
            generate_html_report(report, [f4, f5], ["Prior Contributions", "Posterior Contributions"], resultat, data)

    # If there are no observations:
    else:

        # If the inputs are a specific specie with a specific MeSH, we need export the data table
        if update_data:

            # Rename weight column
            data.rename(columns={'weights':'PriorWeights'}, inplace = True)

            # As null weight have been removed during computation, we use SPECIE instead of index as key
            data["PriorLogOdds"] = np.NaN
            data["PriorLog2FC"] = np.NaN

            for j in data.index:
                ctb_cdf = ss.beta.cdf(p, prior_mix.alpha[j], prior_mix.beta[j])
                data.loc[j, "PriorLogOdds"] = compute_log_odds(ctb_cdf)
                data.loc[j, "PriorLog2FC"] = np.log2(prior_mix.l_mu[j]/p)

        resultat = dict(zip(["TOTAL_PMID_SPECIE", "COOC", "Mean", "CDF", "LogOdds", "Log2FC", "priorLogOdds", "priorLog2FC", "NeighborhoodInformation"], [n, k, prior_mix.mu, prior_mix_CDF, prior_log_odds, prior_log2fc, np.NaN, np.NaN, True]))

        # Plot figure ? Only when the inputs are a specific specie with a specific MeSH
        if report:

            # If names have been provided, use them instead of species labels in Figures:
            if "SPECIE_NAME" in data.columns:
                names = data["SPECIE_NAME"].tolist()
            else:
                names = data["SPECIE"]

            # We set the top to top 10:
            top = 10

            # We select the union of the top 10 contributors in the both groups and then assign a unique color to it in a dict (so max number of contributors displayed is 20)
            set_contributors = set([names[i] for i in np.argsort(prior_mix.weights)[::-1][:top]])

            # We need to keep the same color palette between the both plots
            palette = dict(zip(set_contributors, cm.tab20(np.linspace(0, 1, len(set_contributors)))))

            # Plot
            # f1 = plot_mix_distributions_plotly(prior_mix, names, seq, "Neighbourhood components", palette, top)

            # Contribution plot
            f2 = contributions_plot(data, names, "PriorWeights", "PriorLogOdds")

            # Generate report:
            generate_html_report(report, [f2], ["Prior Contributions"], resultat, data)

    return resultat


def specie2mesh(index, table_cooc, table_species_corpora, weights, table_mesh, forget):
    """This function is used to computed associations from a specific specie against all available MeSHs.

    Args:
        index (int): index of the specie in the metabolic network
        table_cooc (pandas.DataFrame): table of co-occurences
        table_species_corpora (pandas.DataFrame): table of specie corpora
        weights (numpy): weight matrix
        table_mesh (pandas.DataFrame): table of MeSH corpora
        forget (bool): Keep only prior information from the neighborhood, removing specie's observation

    Returns:
        [pd.DataFrame]: association table
    """

    # Create result Dataframe from MeSH list
    mesh_list = table_mesh["MESH"].tolist()
    indexes = range(0, len(mesh_list))
    store = []
    
    # Prepare data table
    table_species_corpora["weights"] = weights[:, index].tolist()

    with progressbar.ProgressBar(max_value=len(indexes)) as bar:
        for i in indexes:
            mesh = mesh_list[i]
            # Get cooc vector. It only contains species that have at least one article, need to left join.
            cooc = table_cooc[table_cooc["MESH"] == mesh][["SPECIE", "COOC"]]
            # Get data
            data = pd.merge(table_species_corpora, cooc, on="SPECIE", how="left").fillna(0)

            # If forget option is true, remove observation of the studied specie
            if forget:
                data.loc[index, ["TOTAL_PMID_SPECIE", "COOC"]] = [0, 0]

            # Get MeSH info
            MeSH_info = table_mesh[table_mesh["MESH"] == mesh]
            p = float(MeSH_info["P"])

            # Computation
            r = computation(index, data, p, float(MeSH_info["alpha_prior"]), float(MeSH_info["beta_prior"]), seq=0.0001)
            store.append(r)
            bar.update(i)

    df_ = pd.DataFrame(store)
    df_.insert(0, "MESH", mesh_list)
    return df_

def mesh2specie(mesh, table_cooc, table_species_corpora, weights, table_mesh, forget):
    """This function is used to computed associations from a specific MeSH against all available species.

    Args:
        index (int): index of the specie in the metabolic network
        table_cooc (pandas.DataFrame): table of co-occurences
        table_species_corpora (pandas.DataFrame): table of specie corpora
        weights (numpy): weight matrix
        table_mesh (pandas.DataFrame): table of MeSH corpora
        forget (bool): Keep only prior information from the neighborhood, removing specie's observation.

    Returns:
        [pd.DataFrame]: association table
    """
    specie_list = table_species_corpora["SPECIE"].tolist()
    indexes = range(0, len(specie_list))
    store = []

    # Get MeSH info
    MeSH_info = table_mesh[table_mesh["MESH"] == mesh]
    cooc = table_cooc[table_cooc["MESH"] == mesh][["SPECIE", "COOC"]]
    p = float(MeSH_info["P"])

    # Browser all species
    with progressbar.ProgressBar(max_value=len(indexes)) as bar:
        for i in indexes:
            table_species_corpora["weights"] = weights[:, i].tolist()
            data = pd.merge(table_species_corpora, cooc, on="SPECIE", how="left").fillna(0)

            # If forget option is true, remove observation of the studied specie
            if forget:
                data.loc[i, ["TOTAL_PMID_SPECIE", "COOC"]] = [0, 0]

            # Computation
            r = computation(i, data, p, float(MeSH_info["alpha_prior"]), float(MeSH_info["beta_prior"]), seq=0.0001)
            store.append(r)
            bar.update(i)

    df_ = pd.DataFrame(store)
    df_.insert(0, "SPECIE", specie_list)
    return df_

def association_file(f, table_cooc, table_species_corpora, weights, table_mesh, forget, species_name_path, out):
    """This function is used to compute all associations specified in a Dataframe (SPECIE, MESH)

    Args:
        f (pandas.Dataframe): A two columns Dataframe storing all SPECIE - MESH pairs that need to be computed.
        table_cooc (pandas.DataFrame): table of co-occurences
        table_species_corpora (pandas.DataFrame): table of specie corpora
        weights (numpy): weight matrix
        table_mesh (pandas.DataFrame): table of MeSH corpora
        forget (bool): Keep only prior information from the neighborhood, removing specie's observation.

    Returns:
        [pd.DataFrame]: association table
    """

    # Test if all provided species and MeSH are present in the network :*
    if (sum(~f["SPECIE"].isin(table_species_corpora["SPECIE"]))) or (sum(~f["MESH"].isin(table_mesh["MESH"]))):
        print("It seems that some species or MeSH descriptors provided in the file are not present in the graph. They will be removed !")
        unavailable = f.loc[~f["SPECIE"].isin(table_species_corpora["SPECIE"]) | ~f["MESH"].isin(table_mesh["MESH"])]
        print("Rows with unavailable info: \n" + unavailable.to_string())
        # Keep only rows with available info:
        f = f.loc[f["SPECIE"].isin(table_species_corpora["SPECIE"]) & f["MESH"].isin(table_mesh["MESH"])]

    store = []
    n = f.shape[0]
    
    # Add names
    if species_name_path:
        table_species_corpora = add_names(table_species_corpora, species_name_path, None)
        table_species_corpora = table_species_corpora[["SPECIE_NAME", "SPECIE", "TOTAL_PMID_SPECIE"]]
    
    # Browse associations
    with progressbar.ProgressBar(max_value=n) as bar:
        for i in range(0, n):
            specie = str(f.iloc[[i], 0].item())
            mesh = str(f.iloc[[i], 1].item())
            index = table_species_corpora[table_species_corpora["SPECIE"] == specie].index[0]
            
            # Prepare data
            table_species_corpora["weights"] = weights[:, index].tolist()
            cooc = table_cooc[table_cooc["MESH"] == mesh][["SPECIE", "COOC"]]
            data = pd.merge(table_species_corpora, cooc, on="SPECIE", how="left").fillna(0)

            # If forget option is true, remove observation of the studied specie
            if forget:
                data.loc[index, ["TOTAL_PMID_SPECIE", "COOC"]] = [0, 0]

            # Get MeSH info
            MeSH_info = table_mesh[table_mesh["MESH"] == mesh]
            p = float(MeSH_info["P"])

            # html report
            path_report = os.path.join(out, "report_" + specie + "_" + mesh + ".html")

            # Computation
            r = computation(index, data, p, float(MeSH_info["alpha_prior"]), float(MeSH_info["beta_prior"]), seq=0.0001, update_data=True, report=path_report, species_name_path=species_name_path)

            out_data = os.path.join(out, "data_" + specie + "_" + mesh + ".csv")
            data.to_csv(out_data, index=False)

            store.append(r)
            bar.update(i)

    associations = pd.concat([f, pd.DataFrame(store)], axis=1)
    return associations

def add_names(result, species_name_path, mesh_name_path):
    #TODO faire par index de colonne pour éviter de devoir spécifier SPECIE_NAME en header du fichier
    # Add species' names if provided
    if species_name_path:
        species_name = pd.read_csv(species_name_path)
        result = pd.merge(result, species_name, on="SPECIE", how="left")
    # Add MESHs' names if provided
    if mesh_name_path:
        mesh_name = pd.read_csv(mesh_name_path)
        result = pd.merge(result, mesh_name, on="MESH", how="left")
    return result
