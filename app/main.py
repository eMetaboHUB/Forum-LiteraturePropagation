import argparse
import sys
import os
from propagation import *

# Get arguments
# python app/main.py --graph="tests/data/test_urea.gml" --specie.corpora="data/species_cid_pmid.csv" --specie.cooc="data/species_cid_mesh_pmid.csv" --mesh.corpora="data/mesh_pmid.csv"
parser = argparse.ArgumentParser()
parser.add_argument("--graph", help="path to the metabolic network compound graph", type = str, required = True, dest = 'g_path')
parser.add_argument("--specie.corpora", help="path to the species corpus size file ", type = str, required = True, dest = 'specie_corpora_path')
parser.add_argument("--specie.cooc", help="path to the species MeSH co-occurences file ", type = str, required = True, dest = 'specie_mesh_path')
parser.add_argument("--mesh", help="The studied MeSH. This argument is incompatible with the 'file' argument. The program will return all association between this MeSH and all species in the metabolic network. If the 'specie' option is also set, only this association specific association will be computed.", type = str, required = False, dest = 'mesh')
parser.add_argument("--specie", help="The studied specie. This argument is incompatible with the 'file' argument. The program will return all association between this specie and all MeSHs. If the 'mesh' option is also set, only this association specific association will be computed.", type = str, required = False, dest = 'specie')
parser.add_argument("--file", help="Path to a file containing pairs of SPECIE and MESH associations to be computed (format csv: SPECIE, MESH). This argument is incompatible with the 'mesh' and 'specie' arguments.", type = str, required = False, dest = 'file')
parser.add_argument("--forget", help="Only the prior from neighborhood will be used in the computation, observations of treated species are set to null. Default = False", type = bool, default = False, required = False, dest = 'forget')
parser.add_argument("--alpha", help="The damping factor (alpha). It could be a single value, or several values to test different parametrizations. All provided alpha values will be tested against all provided sample size values. Default = 0.1", nargs = "*", type = float, default = [0.1], required = False, dest = 'alpha')
parser.add_argument("--sample_size", help="The sample size parameter. It could be a single value, or several values to test different parametrizations. All provided sample size values will be tested against all provided alpha values. Default = 100", nargs = "*", type = int, default = [100], required = False, dest = 'ss')
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
g = import_metabolic_network(args.g_path)

print("> Import species corpora sizes ... ", end = '')
table_species_corpora = import_and_map_indexes(args.specie_corpora_path, g)

# Compute total number of cpd-articles mentions
TOTAL_CPD_MENTIONS = table_species_corpora['TOTAL_PMID_SPECIE'].sum()
print("Ok")

print("> Import species-MeSH co-occurences ... ", end = '')
table_coocurences = import_and_map_indexes(args.specie_mesh_path, g)
print("Ok")

# Compute the total number of mentions between a compound and an article, that also involved MeSHs
table_mesh_corpora = table_coocurences.groupby('MESH', as_index=False)[['COOC']].sum().rename(columns={"COOC": "TOTAL_CPD_MENTION_MESH"})

# Compute MeSH probabilities normalising by the total number of cpd-article mentions
table_mesh_corpora["P"] = table_mesh_corpora["TOTAL_CPD_MENTION_MESH"]/TOTAL_CPD_MENTIONS

# Test if provided MeSH exists
if args.mesh and (not args.mesh in table_mesh_corpora["MESH"].tolist()):
    print("Unknown MeSH identifier: " + args.mesh)
    sys.exit(3)

# Test if provided specie exists
if args.specie and (not args.specie in table_species_corpora["SPECIE"].tolist()):
    print("Unknown specie identifier: " + args.specie)
    sys.exit(3)

# Test if provided file exists
if args.file:
    try:
        f = pd.read_csv(args.file)
    except Exception as e:
        print("Error while trying to read association file. \n" + str(e))


alpha_set = args.alpha
sample_size_set = args.ss

if 0 in sample_size_set:
    print("\n /!\ 0 is not allowed for sample_size.")
    sample_size_set.remove(0)

# Compute analysis

print("\nParameters:\n")
print("- forget: " + str(args.forget))
print("- damping factor (alpha): " + str(alpha_set))
print("- sample size: " + str(sample_size_set))


for alpha in alpha_set:

    # Compute network analysis
    print("\n- Compute weights using alpha = " + str(alpha))
    probabilities = propagation_volume(g, alpha = alpha)
    weights = compute_weights(probabilities, table_species_corpora)
    
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
            f_out_name = os.path.splitext(os.path.basename(args.file))[0]
            out = os.path.join(out_path, f_out_name + "_" + str(alpha) + "_" + str(sample_size) + ("_Forget" * args.forget) + ".csv")
            print("Export results in " + out)
            r.to_csv(out, index = False)

        # If only a specie has been provided
        elif args.specie and not args.mesh:
            print("\nCompute associations between " + args.specie + " and all MeSHs")
            index = int(table_species_corpora[table_species_corpora["SPECIE"] == args.specie]["index"])
            r = specie_mesh(index, table_coocurences, table_species_corpora, weights, table_mesh_corpora_work, args.forget)
            out = os.path.join(out_path, args.specie + "_" + str(alpha) + "_" + str(sample_size) + ("_Forget" * args.forget) + ".csv")
            print("Export results in " + out)
            r.to_csv(out, index = False)
        
        # If only a mesh has been provided
        elif args.mesh and not args.specie:
            print("\nCompute associations between " + args.mesh + " and all species")
            r = mesh_specie(args.mesh, table_coocurences, table_species_corpora, weights, table_mesh_corpora_work, args.forget)
            out = os.path.join(out_path, args.mesh + "_" + str(alpha) + "_" + str(sample_size) + ("_Forget" * args.forget) + ".csv")
            print("Export results in " + out)
            r.to_csv(out, index = False)

        elif args.mesh and args.specie:
            print("\nCompute associations between " + args.specie + "and" + args.mesh)
            index = int(table_species_corpora[table_species_corpora["SPECIE"] == args.specie]["index"])
            table_species_corpora["weights"] = weights[:, index].tolist()
            cooc = table_coocurences[table_coocurences["MESH"] == args.mesh][["index", "COOC"]]
            data = pd.merge(table_species_corpora, cooc, on = "index", how = "left").fillna(0)
            if args.forget:
                data.loc[data["index"] == index, ["TOTAL_PMID_SPECIE", "COOC"]] = [0, 0]
            MeSH_info = table_mesh_corpora_work[table_mesh_corpora_work["MESH"] ==  args.mesh]
            p = float(MeSH_info["P"])
            res = computation(index, data, p, float(MeSH_info["alpha_prior"]), float(MeSH_info["beta_prior"]), seq = 0.0001, plot = True)
            df_ = pd.DataFrame({"SPECIE": args.specie, "MESH": args.mesh, "Mean": [res.Mean], "CDF": [res.CDF], "Log2FC": [res.Log2FC], "priorCDFratio": [res.priorCDFratio], "Score": [res.Score]})
            out = os.path.join(out_path, args.specie + "_" + args.mesh + "_" + str(alpha) + "_" + str(sample_size) + ("_Forget" * args.forget) + ".csv")
            print("Export results in " + out)
            df_.to_csv(out, index = False)
        
        else:
            print("Nothing to do ...")









