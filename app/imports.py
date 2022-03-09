import os
import sys
import igraph as ig
import pandas as pd
import numpy as np



def import_metabolic_network(path, undirected=True, format="gml"):
    """This function is used to import an metabolic network, by default as undirected

    Args:
        path (str): path to the metabolic
        undirected (bool, optional): Graph directed or undirected. Defaults to True.
        format (str, optional): Graph format. Defaults to gml.

    Returns:
        [igraph.Graph]: The compound graph
    """
    if not os.path.isfile(path):
        print("[ERROR] File not found at " + path)
        return None

    try:
        g = ig.read(path, format=format)

    except Exception as e:
        print("[ERROR] While reading graph file at " + path + ": " + str(e))
        return None

    print("[INFO] Metabolic network has been imported successfully.")
    print("[INFO] Number of vertices: " + str(g.vcount()))
    print("[INFO] Number of edges: " + str(g.ecount()))

    if undirected:
        print("[INFO] Used as undirected")
        g.to_undirected()

    else:
        print("[INFO] Used as directed")

    print("[INFO] Extract largest component")
    g = g.clusters().giant()

    # Test if g is a multiple graph, if so sum edges weights:
    if g.has_multiple():
        print("[WARNING] The graph has multiple edges, merge 'Weight' attribute by sum.")
        g.simplify(multiple=True, combine_edges=dict(Weight="sum"))

    # Test if the graph is connected
    if not g.is_connected():
        print("[INFO] The graph needs to be connected: there is a path from any point to any other point in the graph")
        sys.exit(1)

    return g

def import_table(path):
    """This function is a method to read a table. It check if the file exists, and handle exceptions while importing it with pandas.

    Args:
        path (string): path to the table

    Returns:
        (pandas.DataFrame): The imported table
    """
    try:
        data = pd.read_csv(path)

    except pd.errors.ParserError as except_parsing_error:
        print("[ERROR] " + path + " has incorrect format.\n" + str(except_parsing_error))
        return None

    except FileNotFoundError:
        print("[ERROR] File not found at " + path)
        return None

    return data

def import_and_map(path, g, id_att):
    """This function is used to import tables related to species. It also add a new column containing the species index in the graph g.
    Args:
        path (string): path to the table to import
        g (igraph.Graph): The compound graph
        id_att (str): The name of the vertex attribute containing the specie identifier (eg. M_m0001c) in the graph.

    Returns:
        (pandas.DataFrame): The imported tabke, keeping species that are present in the graph
    """

    # Create a table to map species labels (SPECIE column) to indexes in the graph.
    graph_table = pd.DataFrame({"SPECIE": g.vs[id_att]})
    data = import_table(path)

    # If data or graph have not been well imported, return None
    if (data is None) or (g is None):
        return None

    # Test if merge is possible :
    if not list(set(data["SPECIE"]).intersection(graph_table["SPECIE"])):
        print("[ERROR] No 'SPECIE' identifiers in the raw data (" + path + ") is matching with identifiers in the graph (id_att = " + id_att + ").")
        return None

    # Merge
    coocurences = pd.merge(graph_table, data, on="SPECIE", how="left")

    return coocurences
