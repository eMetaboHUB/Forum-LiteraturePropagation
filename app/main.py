import argparse
import sys
import os
from propagation import *
from imports import import_metabolic_network, import_and_map
from weights import create_probabilities_and_weights

# Get arguments
# python app/main.py --graph="data/Recon2_compound-graph.gml" --undirected --specie.corpora="data/species_cid_pmid.csv" --specie.cooc="data/species_cid_mesh_pmid.csv" --specie="M_CE6251" --mesh="D011565" --out="data/tests"  --alpha 0.4 --sample_size 100
parser = argparse.ArgumentParser()
parser.add_argument("--graph", help="Path to the metabolic network compound graph", type=str, required=True, dest='g_path')
parser.add_argument("--undirected", help="Is the graph undirected ? If the graph is directed the probability are extracted from the 'Weight' attribute of the edges. In case of a multi-graph, the multiple edges are merged by summing the weights. If the graph in undirected transition probabilities are computed form the adjacency matrix.", action='store_true')
parser.add_argument("--specie.corpora", help="Path to the SPECIE corpus file ", type=str, required=True, dest='specie_corpora_path')
parser.add_argument("--specie.cooc", help="Path to the SPECIE-MeSH co-occurences file ", type=str, required=True, dest='specie_mesh_path')
parser.add_argument("--mesh", help="The studied MeSH descriptor identifier (eg. D010661). This argument is incompatible with the 'file' argument. This will return all associations between this MeSH and all species in the metabolic network. If the 'specie' option is also set, only this specific association will be computed.", type=str, required=False, dest='mesh')
parser.add_argument("--specie", help="TThe studied specie identifier in the network  according to the *id_att* value (eg. M_m01778c). This argument is incompatible with the 'file' argument. The program will return all associations between this specie and all MeSHs. If the 'mesh' option is also set, only this association specific association will be computed.", type=str, required=False, dest='specie')
parser.add_argument("--file", help="Path to a file containing pairs of SPECIE and MESH associations to be computed (format csv: SPECIE, MESH). This argument is incompatible with the 'mesh' and 'specie' arguments.", type=str, required=False, dest='file')
parser.add_argument("--forget", help="Only the prior from neighborhood will be used in the computation, observations of treated species are 'forgot'.", action='store_true')
parser.add_argument("--alpha", help="The damping factor (alpha). It could be a single value, or several values to test different parametrizations. All provided alpha values will be tested against all provided sample size values. Default = 0.1", nargs="*", type=float, default=[0.1], required=False, dest='alpha')
parser.add_argument("--sample_size", help="The sample size parameter. It could be a single value, or several values to test different parametrizations. All provided sample size values will be tested against all provided alpha values. Default = 100", nargs="*", type=int, default=[100], required=False, dest='ss')
parser.add_argument("--q", help="The tolerance threshold for PPR probabilities. If q is negative the default value is used. The Default = 1/(N - 1)", type=float, default=-1, required=False, dest='q')
parser.add_argument("--id_att", help=" The name of the vertex attribute containing the SPECIE identifier (eg. M_m0001c) in the graph. These identifiers must match with those provided in the 'SPECIE' column of specie.corpora and specie.cooc files. Default is the 'label' attribute", type=str, default="label", required=False, dest='id_att')
parser.add_argument("--species_name_path", help="Path to the file containing species names in a .csv format. First column should be named 'SPECIE' and contains the SPECIE identifiers (eg. M_m0001c) and the second column should be named 'SPECIE_NAME' and contains the species' chemical names", type=str, required=False, dest='species_name_path')
parser.add_argument("--meshs_name_path", help="Path to the file containing species names in a .csv format. First column should be named 'MESH' and contains the MESH identifiers (eg. D000017) and the second column should be named 'MESH_NAME' and contains the MeSH labels", type=str, required=False, dest='meshs_name_path')
parser.add_argument("--out", help="Path to the output directory", type=str, required=True, dest='out')

args = parser.parse_args()

# Parse options
if args.file and (args.specie or args.mesh):
    print("[INFO] The 'file' option and 'mesh' and/or 'specie' are incomptibles: \n")
    print("(1) To compute associations between a specific specie and all MeSH, use only the 'specie' option")
    print("(2) To compute associations between a specific MeSH and all specie, use only the 'mesh' option")
    print("(3) To compute a specific association, set both the 'specie' and the 'mesh' option")
    print("(4) To compute a specific set of associations provided as pairs in a file, use the file option")
    sys.exit(2)

out_path = args.out

# Import data
g_name = os.path.splitext(os.path.basename(args.g_path))[0]
g = import_metabolic_network(args.g_path, undirected=args.undirected)
if g is None:
    print("[INFO] Exit due to errors during data import")
    sys.exit(1)

print("[INFO] Import species corpora sizes ... ", end='')
table_species_corpora = import_and_map(args.specie_corpora_path, g, args.id_att)
if table_species_corpora is None:
    print("[INFO] Exit due to errors during data import")
    sys.exit(1)

# Fill na if some species are not indicated in the file
table_species_corpora = table_species_corpora.fillna(0)

# Compute total number of cpd-articles mentions
total_cpd_mention = table_species_corpora['TOTAL_PMID_SPECIE'].sum()
print("Ok")

print("[INFO] Import species-MeSH co-occurences ... ", end='')
table_coocurences = import_and_map(args.specie_mesh_path, g, args.id_att)
if table_coocurences is None:
    print("[INFO] Exit due to errors during data import")
    sys.exit(1)
print("Ok")

# Compute the total number of mentions between a compound and an article, that also involved MeSHs
table_mesh_corpora = table_coocurences.groupby('MESH', as_index=False)[['COOC']].sum().rename(columns={"COOC": "TOTAL_CPD_MENTION_MESH"})

# Compute MeSH probabilities normalising by the total number of cpd-article mentions
table_mesh_corpora["P"] = table_mesh_corpora["TOTAL_CPD_MENTION_MESH"]/total_cpd_mention

# Test if provided MeSH exists
if args.mesh and (not args.mesh in table_mesh_corpora["MESH"].tolist()):
    print("[INFO] Unknown MeSH identifier: " + args.mesh)
    sys.exit(1)

# Test if provided specie exists
if args.specie and (not args.specie in table_species_corpora["SPECIE"].tolist()):
    print("[INFO] Unknown specie identifier: " + args.specie)
    sys.exit(1)

# Test if provided file exists
if args.file:

    try:
        f = pd.read_csv(args.file)

    except pd.errors.ParserError as except_parsing_error:
        print("[ERROR] " + args.file + " has incorrect format.\n" + str(except_parsing_error))
        sys.exit(1)

    except FileNotFoundError:
        print("[ERROR] File not found at " + args.file + "\n")
        sys.exit(1)

    # Test columns:
    if (set(f.columns.to_list()) != set(["SPECIE", "MESH"])) or f.empty:
        print("[ERROR] Bad formating for association file. The file need to contain only 2 columns: SPECIE and MESH")
        sys.exit(1)



alpha_set = args.alpha
sample_size_set = args.ss
q = args.q

if q <= 0:
    print("[INFO] Set tolerance threshold to default value")
    q = 1/(len(g.vs) - 1)

if 0 in sample_size_set:
    print("[WARNING] 0 is not allowed for sample_size.")
    sample_size_set.remove(0)

# Compute analysis

print("[INFO] Parameters:\n")
print("- forget: " + str(args.forget))
print("- damping factor (alpha): " + str(alpha_set))
print("- sample size: " + str(sample_size_set))
print("- q: " + str(q))

# Extract labels for supplementary tables
l = table_species_corpora["SPECIE"]

for alpha in alpha_set:

    # Compute network analysis
    probabilities, weights = create_probabilities_and_weights(g, g_name, alpha, table_species_corpora, q, args.id_att)
    df_Entropy = compute_entropy_matrix(weights, l)
    df_contributors_distances = compute_contributors_distances(weights, g, l)
    df_contributors_corpora_sizes = compute_contributors_corpora_sizes(weights, table_species_corpora, l)
    df_nb_ctbs = compute_contributors_number(weights, l)

    for sample_size in sample_size_set:
        print("[INFO] Compute MeSH priors using sample size = " + str(sample_size))
        # Compute prior parameters:
        mesh_priors = table_mesh_corpora["P"].apply(estimate_prior_distribution_mesh, sample_size=sample_size)
        mesh_priors = pd.DataFrame(mesh_priors.tolist(), columns=["alpha_prior", "beta_prior"])
        table_mesh_corpora_work = pd.concat([table_mesh_corpora, mesh_priors], axis=1)

        print("[INFO] Treating alpha = " + str(alpha) + " and sample_size = " + str(sample_size))

        # If an SPECIE-MESH file was provided:

        if args.file and (not f.empty):

            # Compute predictions
            print("[INFO] Compute association from file: " + args.file)
            r = association_file(f, table_coocurences, table_species_corpora, weights, table_mesh_corpora_work, args.forget, args.species_name_path, args.out)

            # Add diagnostic values for each prior
            r = pd.merge(r, df_Entropy, on="SPECIE", how="left")
            r = pd.merge(r, df_contributors_distances, on="SPECIE", how="left")
            r = pd.merge(r, df_contributors_corpora_sizes, on="SPECIE", how="left")
            r = pd.merge(r, df_nb_ctbs, on="SPECIE", how="left")

            # Export
            f_out_name = os.path.splitext(os.path.basename(args.file))[0]
            out = os.path.join(out_path, f_out_name + "_" + str(alpha) + "_" + str(sample_size) + ("_Forget" * args.forget) + ".csv")
            r = add_names(r, args.species_name_path, args.meshs_name_path)
            print("[INFO] Export results in " + out)
            r.to_csv(out, index=False)

        # If only a specie has been provided
        elif args.specie and not args.mesh:

            print("[INFO] Compute associations between " + args.specie + " and all MeSHs")

            # Prepare data table
            index = table_species_corpora[table_species_corpora["SPECIE"] == args.specie].index[0]

            # Computation
            r = specie2mesh(index, table_coocurences, table_species_corpora, weights, table_mesh_corpora_work, args.forget)

            # Add diagnostic values for each prior
            r["Entropy"] = float(df_Entropy.loc[index, "Entropy"])
            r["CtbAvgDistance"] = float(df_contributors_distances.loc[index, "CtbAvgDistance"])
            r["CtbAvgCorporaSize"] = float(df_contributors_corpora_sizes.loc[index, "CtbAvgCorporaSize"])
            r["NbCtb"] = float(df_nb_ctbs.loc[index, "NbCtb"])
            out = os.path.join(out_path, args.specie + "_" + str(alpha) + "_" + str(sample_size) + ("_Forget" * args.forget) + ".csv")

            # Export
            r = add_names(r, None, args.meshs_name_path)
            print("[INFO] Export results in " + out)
            r.to_csv(out, index=False)

        # If only a mesh has been provided
        elif args.mesh and not args.specie:

            print("[INFO] Compute associations between " + args.mesh + " and all species")

            # Computation
            r = mesh2specie(args.mesh, table_coocurences, table_species_corpora, weights, table_mesh_corpora_work, args.forget)

            # Add diagnostic values for each prior
            r["Entropy"] = df_Entropy["Entropy"]
            r["CtbAvgDistance"] = df_contributors_distances["CtbAvgDistance"]
            r["CtbAvgCorporaSize"] = df_contributors_corpora_sizes["CtbAvgCorporaSize"]
            r["NbCtb"] = df_nb_ctbs["NbCtb"]

            # Export
            out = os.path.join(out_path, args.mesh + "_" + str(alpha) + "_" + str(sample_size) + ("_Forget" * args.forget) + ".csv")
            r = add_names(r, args.species_name_path, None)
            print("[INFO] Export results in " + out)
            r.to_csv(out, index=False)

        elif args.mesh and args.specie:

            print("[INFO] Compute associations between " + args.specie + " and " + args.mesh)

            # Prepare data table
            index = table_species_corpora[table_species_corpora["SPECIE"] == args.specie].index[0]
            table_species_corpora["weights"] = weights[:, index].tolist()
            cooc = table_coocurences[table_coocurences["MESH"] == args.mesh][["SPECIE", "COOC"]]
            data = pd.merge(table_species_corpora, cooc, on="SPECIE", how="left").fillna(0)

            # Forget specie's literature ?
            if args.forget:
                data.loc[index, ["TOTAL_PMID_SPECIE", "COOC"]] = [0, 0]

            # MeSH proba
            MeSH_info = table_mesh_corpora_work[table_mesh_corpora_work["MESH"] == args.mesh]
            p = float(MeSH_info["P"])
            out_path_report = os.path.join(out_path, "report_" + args.specie + "_" + args.mesh + "_" + str(alpha) + "_" + str(sample_size) + ("_Forget" * args.forget) + ".html")
            
            # Computation
            res = computation(index, data, p, float(MeSH_info["alpha_prior"]), float(MeSH_info["beta_prior"]), seq=0.0001, species_name_path=args.species_name_path, update_data=True, report=out_path_report)
            df_ = pd.concat([pd.DataFrame([{"SPECIE": args.specie, "MESH": args.mesh}]), pd.DataFrame([res])], axis=1)

            # Add diagnostic values for each prior
            df_["Entropy"] = df_Entropy.loc[index, "Entropy"]
            df_["CtbAvgDistance"] = df_contributors_distances.loc[index, "CtbAvgDistance"]
            df_["CtbAvgCorporaSize"] = df_contributors_corpora_sizes.loc[index, "CtbAvgCorporaSize"]
            df_["NbCtb"] = df_nb_ctbs.loc[index, "NbCtb"]
            out = os.path.join(out_path, args.specie + "_" + args.mesh + "_" + str(alpha) + "_" + str(sample_size) + ("_Forget" * args.forget) + ".csv")

            # Export
            df_ = add_names(df_, args.species_name_path, args.meshs_name_path)
            data = add_names(data, args.species_name_path, None)
            print("[INFO] Export results in " + out)
            df_.to_csv(out, index=False)
            out_data = os.path.join(out_path, "data_" + args.specie + "_" + args.mesh + "_" + str(alpha) + "_" + str(sample_size) + ("_Forget" * args.forget) + ".csv")
            data.to_csv(out_data, index=False)
        else:
            print("Nothing to do ...")
