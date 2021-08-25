import argparse
import sys
import os
from propagation import *

# Get arguments
# python app/main.py --graph="data/Recon2_compound-graph.gml" --undirected --specie.corpora="data/species_cid_pmid.csv" --specie.cooc="data/species_cid_mesh_pmid.csv" --specie="M_CE6251" --mesh="D011565" --out="data/tests"  --alpha 0.4 --sample_size 100
parser = argparse.ArgumentParser()
parser.add_argument("--graph", help="path to the metabolic network compound graph", type = str, required = True, dest = 'g_path')
parser.add_argument("--undirected", help="Is the graph undirected ?(Cf. doc)", action='store_true')
parser.add_argument("--specie.corpora", help="path to the species corpus size file ", type = str, required = True, dest = 'specie_corpora_path')
parser.add_argument("--specie.cooc", help="path to the species MeSH co-occurences file ", type = str, required = True, dest = 'specie_mesh_path')
parser.add_argument("--mesh", help="The studied MeSH. This argument is incompatible with the 'file' argument. The program will return all association between this MeSH and all species in the metabolic network. If the 'specie' option is also set, only this association specific association will be computed.", type = str, required = False, dest = 'mesh')
parser.add_argument("--specie", help="The studied specie. This argument is incompatible with the 'file' argument. The program will return all association between this specie and all MeSHs. If the 'mesh' option is also set, only this association specific association will be computed.", type = str, required = False, dest = 'specie')
parser.add_argument("--file", help="Path to a file containing pairs of SPECIE and MESH associations to be computed (format csv: SPECIE, MESH). This argument is incompatible with the 'mesh' and 'specie' arguments.", type = str, required = False, dest = 'file')
parser.add_argument("--forget", help="Only the prior from neighborhood will be used in the computation, observations of treated species are 'forgot'.", action='store_true')
parser.add_argument("--alpha", help="The damping factor (alpha). It could be a single value, or several values to test different parametrizations. All provided alpha values will be tested against all provided sample size values. Default = 0.1", nargs = "*", type = float, default = [0.1], required = False, dest = 'alpha')
parser.add_argument("--sample_size", help="The sample size parameter. It could be a single value, or several values to test different parametrizations. All provided sample size values will be tested against all provided alpha values. Default = 100", nargs = "*", type = int, default = [100], required = False, dest = 'ss')
parser.add_argument("--q", help="The tolerance threshold for PPR probabilities. If q is negative the default value is used. The Default = 1/(N  -1)", type = float, default = -1, required = False, dest = 'q')
parser.add_argument("--id_att", help=" The name of the vertex attribute containing the SPECIE identifier (eg. M_m0001c) in the graph. These identifiers must match with those provided in the 'SPECIE' column of specie.corpora and specie.cooc files. Default is the 'label' attribute", type = str, default = "label", required = False, dest = 'id_att')
parser.add_argument("--species_name_path", help="path to the file containing species names in a .csv format. First column should be named 'SPECIE' and contains the SPECIE identifiers (eg. M_m0001c) and the second column should be named 'SPECIE_NAME' and contains the species' chemical names", type = str, required = False, dest = 'species_name_path')
parser.add_argument("--meshs_name_path", help="path to the file containing species names in a .csv format. First column should be named 'MESH' and contains the MESH identifiers (eg. D000017) and the second column should be named 'MESH_NAME' and contains the MeSH labels", type = str, required = False, dest = 'meshs_name_path')
parser.add_argument("--out", help="path to the output directory", type = str, required = True, dest = 'out')


args = parser.parse_args()

# Parse options
if args.file and (args.specie or args.mesh):
    print("\nThe 'file' option and 'mesh' and/or 'specie' are incomptibles: \n")
    print("(1) To compute associations between a specific specie and all MeSH, use only the 'specie' option")
    print("(2) To compute associations between a specific MeSH and all specie, use only the 'mesh' option")
    print("(3) To compute a specific association, set both the 'specie' and the 'mesh' option")
    print("(4) To compute a specific set of associations provided as pairs in a file, use the file option")
    sys.exit(2)

out_path = args.out

# Import data
g_name = os.path.splitext(os.path.basename(args.g_path))[0]
g = import_metabolic_network(args.g_path, undirected = args.undirected)
if g is None:
    print("\n /!\ Exit due to errors during data import")
    sys.exit(1)

print("> Import species corpora sizes ... ", end = '')
table_species_corpora = import_and_map_indexes(args.specie_corpora_path, g, args.id_att)
if table_species_corpora is None:
    print("\n /!\ Exit due to errors during data import")
    sys.exit(1)

# Fill na if some species are not indicated in the file
table_species_corpora = table_species_corpora.fillna(0)

# Compute total number of cpd-articles mentions
TOTAL_CPD_MENTIONS = table_species_corpora['TOTAL_PMID_SPECIE'].sum()
print("Ok")

print("> Import species-MeSH co-occurences ... ", end = '')
table_coocurences = import_and_map_indexes(args.specie_mesh_path, g, args.id_att)
if table_coocurences is None:
    print("\n /!\ Exit due to errors during data import")
    sys.exit(1)
print("Ok")

# Compute the total number of mentions between a compound and an article, that also involved MeSHs
table_mesh_corpora = table_coocurences.groupby('MESH', as_index=False)[['COOC']].sum().rename(columns={"COOC": "TOTAL_CPD_MENTION_MESH"})

# Compute MeSH probabilities normalising by the total number of cpd-article mentions
table_mesh_corpora["P"] = table_mesh_corpora["TOTAL_CPD_MENTION_MESH"]/TOTAL_CPD_MENTIONS

# Test if provided MeSH exists
if args.mesh and (not args.mesh in table_mesh_corpora["MESH"].tolist()):
    print("Unknown MeSH identifier: " + args.mesh)
    sys.exit(1)

# Test if provided specie exists
if args.specie and (not args.specie in table_species_corpora["SPECIE"].tolist()):
    print("Unknown specie identifier: " + args.specie)
    sys.exit(1)

# Test if provided file exists
if args.file:
    try:
        f = pd.read_csv(args.file)
    except Exception as e:
        print("Error while trying to read association file. \n" + str(e))
        sys.exit(1)
    # Test columns:
    if (set(f.columns.to_list()) != set(["SPECIE", "MESH"])) or f.empty:
        print("Bad formating for association file. The file need to contain only 2 column: SPECIE and MESH")
        sys.exit(1)



alpha_set = args.alpha
sample_size_set = args.ss
q = args.q

if q <= 0:
    print("\n Set tolerance threshold to default value")
    q = 1/(len(g.vs) - 1)

if 0 in sample_size_set:
    print("\n /!\ 0 is not allowed for sample_size.")
    sample_size_set.remove(0)

# Compute analysis

print("\nParameters:\n")
print("- forget: " + str(args.forget))
print("- damping factor (alpha): " + str(alpha_set))
print("- sample size: " + str(sample_size_set))
print("- q: " + str(q))

# Extract labels for supplementary tables
l = table_species_corpora["SPECIE"]

for alpha in alpha_set:

    # Compute network analysis
    probabilities, weights = create_probabilities_and_weights(g, g_name, alpha, table_species_corpora, q, args.id_att)
    df_Entropy = compute_Entropy_matrix(weights, l)
    df_contributors_distances = compute_contributors_distances(weights, g, l)
    df_contributors_corpora_sizes = compute_contributors_corpora_sizes(weights, table_species_corpora, l)
    df_nb_ctbs = compute_contributors_number(weights, l)
    
    for sample_size in sample_size_set:
        print("\n- Compute MeSH priors using sample size = " + str(sample_size))
        # Compute prior parameters:
        mesh_priors = table_mesh_corpora["P"].apply(estimate_prior_distribution_mesh_V2, sample_size = sample_size)
        mesh_priors = pd.DataFrame(mesh_priors.tolist(), columns = ["alpha_prior", "beta_prior"])
        table_mesh_corpora_work = pd.concat([table_mesh_corpora, mesh_priors], axis = 1)
        
        print("\nTreating alpha = " + str(alpha) + " and sample_size = " + str(sample_size))

        # If an SPECIE-MESH file was provided:
        if args.file and (not f.empty):
            print("\nCompute association from file: " + args.file)
            r = association_file(f, table_coocurences, table_species_corpora, weights, table_mesh_corpora_work, args.forget)
            # left join to add entropy, CtbAvgDistance, CtbAvgCorporaSize
            r = pd.merge(r, df_Entropy, on = "SPECIE", how = "left")
            r = pd.merge(r, df_contributors_distances, on = "SPECIE", how = "left")
            r = pd.merge(r, df_contributors_corpora_sizes, on = "SPECIE", how = "left")
            r = pd.merge(r, df_nb_ctbs, on = "SPECIE", how = "left")
            f_out_name = os.path.splitext(os.path.basename(args.file))[0]
            out = os.path.join(out_path, f_out_name + "_" + str(alpha) + "_" + str(sample_size) + ("_Forget" * args.forget) + ".csv")
            r = add_names(r, args.species_name_path, args.meshs_name_path)
            print("Export results in " + out)
            r.to_csv(out, index = False)

        # If only a specie has been provided
        elif args.specie and not args.mesh:
            print("\nCompute associations between " + args.specie + " and all MeSHs")
            index = int(table_species_corpora[table_species_corpora["SPECIE"] == args.specie]["index"])
            r = specie_mesh(index, table_coocurences, table_species_corpora, weights, table_mesh_corpora_work, args.forget)
            # Add Entropy, CtbAvgDistance, CtbAvgCorporaSize, NbCtb of the targeted compound
            r["Entropy"] = float(df_Entropy[df_Entropy["SPECIE"] == args.specie]["Entropy"])
            r["CtbAvgDistance"] = float(df_contributors_distances[df_contributors_distances["SPECIE"] == args.specie]["CtbAvgDistance"])
            r["CtbAvgCorporaSize"] = float(df_contributors_corpora_sizes[df_contributors_corpora_sizes["SPECIE"] == args.specie]["CtbAvgCorporaSize"])
            r["NbCtb"] = float(df_nb_ctbs[df_nb_ctbs["SPECIE"] == args.specie]["NbCtb"])
            out = os.path.join(out_path, args.specie + "_" + str(alpha) + "_" + str(sample_size) + ("_Forget" * args.forget) + ".csv")
            # Add labels to result if provided
            r = add_names(r, None, args.meshs_name_path)
            print("Export results in " + out)
            r.to_csv(out, index = False)
        
        # If only a mesh has been provided
        elif args.mesh and not args.specie:
            print("\nCompute associations between " + args.mesh + " and all species")
            r = mesh_specie(args.mesh, table_coocurences, table_species_corpora, weights, table_mesh_corpora_work, args.forget)
            # Add Entropy, CtbAvgDistance, CtbAvgCorporaSize, NbCtb
            r["Entropy"] = df_Entropy["Entropy"]
            r["CtbAvgDistance"] = df_contributors_distances["CtbAvgDistance"]
            r["CtbAvgCorporaSize"] = df_contributors_corpora_sizes["CtbAvgCorporaSize"]
            r["NbCtb"] = df_nb_ctbs["NbCtb"]
            out = os.path.join(out_path, args.mesh + "_" + str(alpha) + "_" + str(sample_size) + ("_Forget" * args.forget) + ".csv")
            # Add labels to result if provided
            r = add_names(r, args.species_name_path, None)
            print("Export results in " + out)
            r.to_csv(out, index = False)

        elif args.mesh and args.specie:
            print("\nCompute associations between " + args.specie + " and " + args.mesh)
            index = int(table_species_corpora[table_species_corpora["SPECIE"] == args.specie]["index"])
            table_species_corpora["weights"] = weights[:, index].tolist()
            cooc = table_coocurences[table_coocurences["MESH"] == args.mesh][["index", "COOC"]]
            data = pd.merge(table_species_corpora, cooc, on = "index", how = "left").fillna(0)
            if args.forget:
                data.loc[data["index"] == index, ["TOTAL_PMID_SPECIE", "COOC"]] = [0, 0]
            MeSH_info = table_mesh_corpora_work[table_mesh_corpora_work["MESH"] ==  args.mesh]
            p = float(MeSH_info["P"])
            print("P = " + str(p))
            res = computation(index, data, p, float(MeSH_info["alpha_prior"]), float(MeSH_info["beta_prior"]), seq = 0.0001, plot = True, species_name_path = args.species_name_path)
            df_ = pd.DataFrame({"SPECIE": args.specie, "MESH": args.mesh, "TOTAL_PMID_SPECIE": [res.TOTAL_PMID_SPECIE], "COOC": [res.COOC], "Mean": [res.Mean], "CDF": [res.CDF], "Log2FC": [res.Log2FC], "priorCDF": [res.priorCDF], "priorLog2FC": [res.priorLog2FC], "NeighborhoodInformation": [res.NeighborhoodInformation], "Entropy": df_Entropy[df_Entropy["SPECIE"] == args.specie]["Entropy"], "CtbAvgDistance": df_contributors_distances[df_contributors_distances["SPECIE"] == args.specie]["CtbAvgDistance"], "CtbAvgCorporaSize": df_contributors_corpora_sizes[df_contributors_corpora_sizes["SPECIE"] == args.specie]["CtbAvgCorporaSize"], "NbCtb": df_nb_ctbs[df_nb_ctbs["SPECIE"] == args.specie]["NbCtb"]})
            out = os.path.join(out_path, args.specie + "_" + args.mesh + "_" + str(alpha) + "_" + str(sample_size) + ("_Forget" * args.forget) + ".csv")
            # Add labels to result if provided
            df_ = add_names(df_, args.species_name_path, args.meshs_name_path)
            data = add_names(data, args.species_name_path, None)
            print("Export results in " + out)
            df_.to_csv(out, index = False)
            # Export full data
            out_data = os.path.join(out_path, "data_" + args.specie + "_" + args.mesh + "_" + str(alpha) + "_" + str(sample_size) + ("_Forget" * args.forget) + ".csv")
            data.to_csv(out_data, index = False)
            # Vizu:
            # vizu = create_vizu_data(index, probabilities, q, weights, data)
            # out_vizu = os.path.join(out_path, "vizu_" + args.specie + "_" + args.mesh + "_" + str(alpha) + "_" + str(sample_size) + ("_Forget" * args.forget) + ".csv")
            # vizu.to_csv(out_vizu, index = False, sep = "\t")
        else:
            print("Nothing to do ...")









