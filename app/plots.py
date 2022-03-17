import numpy as np
import plotly
import plotly.graph_objects as go
import scipy.stats as ss
import plotly.express as px

from matplotlib.colors import rgb2hex

def plot_distributions_plotly(prior_mix, posterior_mix):
    """This function is used to plot prior distribution against a posterior distribution
    The figure is computed using plotly.

    Args:
        prior_mix (collection): A collection containing information about a prior distribution with x probabilties and associated densities from create_prior_beta_mix or simple_prior with the samping = True
        posterior_mix (collection): A collection containing information about a posterior distribution with x probabilties and associated densities from create_posterior_beta_mix or simple_posterior with the samping = True
    """
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=prior_mix.x,
            y=prior_mix.f,
            line=dict(color="blue"),
            name="Prior"
            )
        )
    fig.add_trace(
        go.Scatter(
            x=posterior_mix.x,
            y=posterior_mix.f,
            line=dict(color="red"),
            name="Posterior"
        )
    )
    fig.update_layout(title="Differences between prior mix and posterior mix distribution",
        xaxis=dict(title="Probability", titlefont_size=25, tickfont_size=20),
        yaxis=dict(title="Density", titlefont_size=25, tickfont_size=20),
        template="simple_white")

    return fig

def plot_mix_distributions_plotly(mix, labels, seq, name, color_palette, top):
    """This function is used to plot distribution of each components of a prior mixture distribution. The function compute itself the densities of each component. Since there could be dozens of contributors, we only plot the top n (top argument) for clarity.
    The figure is computed using plotly.

    Args:
        mix (collections): A collection containing information about the mixture distribution: weights, alpha and beta parameters
        labels (list): A list of compound (or specie) labels associated to each component of the mixture. 
        seq (float): The step used in np.arange to create the x vector of probabilities.
        color_palette (dict):  A dict containing as key the label of the specie and as value a np.array of 4 elements representing the associated color in RGBA format
        top (int): The top n (maximum) of contributors that should be plotted
    """
    fig = go.Figure()
    x = np.arange(0, 1 + seq, seq).tolist()
    weights = mix.weights

    # To plot only the top n contributors, we first order the index of weights in decreasing order
    ordered_i_w = np.flip(np.argsort(weights))

    # We go trough the list of weight until we reach the n'th contributor, or the last contributor if there are les than n
    for i in range(0, min(len(weights), top)):
        it = ordered_i_w[i]
        f = ss.beta.pdf(x, a=mix.alpha[it], b=mix.beta[it])
        y = weights[it] * f
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                line=dict(color=rgb2hex(color_palette[labels[it]])),
                name=labels[it] + ": Beta(" + str(round(mix.alpha[it], 2)) + ", " + str(round(mix.beta[it], 2)) + ") - w = " + str(round(weights[it],2)),
                hovertemplate="<b>Contributor:</b>"+ labels[it] + "<extra></extra>")
            )
    fig.update_layout(title="Top " + str(min(len(weights), top)) + " - " + name + " decomposition",
        xaxis=dict(title="Probability", titlefont_size=25, tickfont_size=20),
        yaxis=dict(title="Density", titlefont_size=25, tickfont_size=20),
        template="simple_white")
    
    return fig


def contributions_plot(data, names, attr, odds_var, limit=0.99):
    """This function is used to produce the contribution plot

    Args:
        data (pandas.DataFrame): Data from computation
        names (list): A list of names
        attr (str): The attribute used to sort contributors (weights or posterioir_weights)
        limit (float, optional): The limit from which contributors will be assigned to the 'others' category by computing the cumulative sum of their posterior weights. Defaults to 0.99.

    Returns:
        [plotly]: The figure
    """
    # First make a copy of data to keep the original clean
    _data = data.copy()
    _data["y"] = "contributors"
    _data["name"] = names
    _data = _data.rename(columns={attr: "Contributions"})
    _data = _data[["SPECIE", "TOTAL_PMID_SPECIE", "COOC", "Contributions", odds_var, "y", "name"]]

    # Sort by posterior weights and determine contributors that belong the the 'others' category
    _data.sort_values(by='Contributions', inplace=True, ascending=False, ignore_index=True)
    cumsum = np.cumsum(_data["Contributions"].tolist())
    l = np.argmin(abs(cumsum - limit))

    # Compute stats for the 'others' category
    others_median = np.median(_data.loc[_data.index[(l + 1):], odds_var])
    others_cooc = np.median(_data.loc[_data.index[(l + 1):], "COOC"])
    others_total_pmid_specie = np.median(_data.loc[_data.index[(l + 1):], "TOTAL_PMID_SPECIE"])

    # Replace 'others' contributors by the 'others' line
    _data.drop(index=_data.index[(l + 1):], inplace=True)
    _data.loc[-1] = ["others", others_total_pmid_specie, others_cooc, (1 - cumsum[l]), others_median, "contributors", "others"]

    # To manage the color scale, we have to make a copy of Log_Odds. As many contributors could have very high or small LogOdds (eg. Inf or -Inf) we restrict their values to a range of -100 - 100 (in Odds, not LogOdds)
    # For contributors that have an Odds higher than 100 or lower than -100, we replace their value by le limit (100 or -100) to restrict the color scale.
    _data["w_LogOdds"] = _data[odds_var]
    _data.loc[_data.w_LogOdds >= np.log(100), "w_LogOdds"] = np.log(100)
    _data.loc[_data.w_LogOdds <= np.log(0.01), "w_LogOdds"] = np.log(0.01)
    _data[odds_var] = [str(v) for v in np.round(_data[odds_var], 2)]

    # hover_dict = dict({"TOTAL_PMID_SPECIE": ":.", "COOC": ":.", "Contributions": ":.2f", "LogOdds": True, "w_LogOdds": False, "y": False})
    hover_dict = dict(zip(("TOTAL_PMID_SPECIE", "COOC", "Contributions", odds_var, "w_LogOdds", "y"), (":.", ":.", ":.2f", True, False, False)))

    fig = px.bar(_data, y="y",
        x="Contributions",
        color="w_LogOdds",
        orientation="h",
        hover_data=hover_dict,
        hover_name="name",
        height=800,
        range_color=[np.log(0.01), np.log(100)],
        facet_col_spacing=1,
        color_continuous_scale=[(0, "blue"), (0.5, "white"), (1, "red")],
        template="seaborn",
        labels={"y": '', "Contributions": "Contributions"})
    
    # The 'len' attribute in important
    fig.update_layout(coloraxis_colorbar=dict(
        title=dict(text="Contributor Odds (in log scale)", font=dict(size=25)),
        tickvals=[np.log(0.01), np.log(0.02), np.log(0.1), 0, np.log(10), np.log(50), np.log(100)],
        ticktext=["<= 0.01", "0.02", "0.1", "0", "10", "50", ">= 100"],
        len=5,
        tickfont=dict(size=20)),
        xaxis=dict(titlefont_size=25, tickfont_size=20),
        yaxis=dict(titlefont_size=25, tickfont_size=20)
    )
    fig.update_traces(marker_line_color='rgb(0,0,0)', marker_line_width=1, opacity=1)
    fig.update_xaxes(range=[0, 1])

    return fig

def generate_html_report(out_path, figs, section_titles, resultat):
    """This function is used to produce the HTML report

    Args:
        out_path (str): the path to write the html report
        figs (list): A list of figures
        section_titles (list): A list of section titles for the figures
        resultat (dict): The result dict with CDF, LogOdds, etc ...
    """
    
    contributors = "<p>".join(["<h2>" + section_titles[i] + "</h2>" + plotly.offline.plot(figs[i], include_plotlyjs=False, output_type='div') for i in range(len(figs))])
    res = "<p>".join(["<b>" + k + ": </b> " + "{:.2e}".format(v) if isinstance(v, float) else "<b>" + k + ": </b> " + str(v) for k, v in resultat.items()])
    html_template = f"""
    <html>
    <head>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    </head>
    <body>
        <h1>NEIGHBORHOOD CONTRIBUTION: HTML REPORT</h1>
        <h2> RESULTS</h2>
        {res}
        <h2> CONTRIBUTORS </h2>
        {contributors}
    </body>
    </html>
    """

    with open(out_path, "w") as f:
        f.write(html_template)
