import argparse
from propagation import *

# Get arguments
# python app/main.py --graph="tests/data/test_urea.gml" --specie.corpora="data/species_cid_pmid.csv" --specie.cooc="data/species_cid_mesh_pmid.csv" --mesh.corpora="data/mesh_pmid.csv"
parser = argparse.ArgumentParser()
parser.add_argument("--graph", help="path to metabolic network compound graph", type = str, required = True, dest = 'g_path')
parser.add_argument("--specie.corpora", help="path to the species corpus size file ", type = str, required = True, dest = 'specie_corpora_path')
parser.add_argument("--specie.cooc", help="path to the species MeSH co-occurences file ", type = str, required = True, dest = 'specie_mesh_path')
parser.add_argument("--mesh.corpora", help="path to the MeSH corpus size file ", type = str, required = True, dest = 'mesh_corpora_path')

args = parser.parse_args()

sample_size = 100
alpha = 0.1

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

# Compute prior parameters:
mesh_priors = table_mesh_corpora["P"].apply(estimate_prior_distribution_mesh_V2, sample_size = sample_size)
mesh_priors = pd.DataFrame(mesh_priors.tolist(), columns = ["alpha_prior", "beta_prior"])

table_mesh_corpora = pd.concat([table_mesh_corpora, mesh_priors], axis = 1)
print("Ok")


mesh = "D052536" # "D002386" # "D018312"
specie = "M_spc_hs" # "M_zymstnl" # "M_tststerone"


# probabilities = propagation_volume(g, alpha = alpha)

# weights = compute_weights(probabilities, table_species_corpora)

if True:
    validation_set = pd.read_csv("data/Validation/validation_set_classic.csv")
    # add results columns
    validation_set = pd.concat([validation_set, pd.DataFrame(columns = ["Mean", "CDF", "Log2FC", "priorCDFratio"])])
    # Compute the total number of mentions between a compound and an article, that also involved MeSHs
    M = table_coocurences.groupby('MESH', as_index=False)[['COOC']].sum().rename(columns={"COOC": "TOTAL_CPD_MENTION_MESH"})
    # Compute MeSH probabilities normalising by the total number of cpd-article mentions
    M["P"] = M["TOTAL_CPD_MENTION_MESH"]/TOTAL_CPD_MENTIONS
    # Iter over associations
    for alpha in [0, 0.1, 0.2, 0.3, 0.4, 0.5]:
        probabilities = propagation_volume(g, alpha = alpha)
        weights = compute_weights(probabilities, table_species_corpora)
        for sample_size in [100000, 1000000]: # 1, 10, 100, 1000, 10000
            # Prepare priors
            mesh_priors = M["P"].apply(estimate_prior_distribution_mesh_V2, sample_size = sample_size)
            mesh_priors = pd.DataFrame(mesh_priors.tolist(), columns = ["alpha_prior", "beta_prior"])
            table_mesh_corpora = pd.concat([M, mesh_priors], axis = 1)
            # Compute validation for set parameters
            print("Treating alpha = " + str(alpha) + " and sample_size = " + str(sample_size))
            for i in range(0, len(validation_set.index)):
                # get row info
                specie = str(validation_set.iloc[[i], 0].item())
                mesh = str(validation_set.iloc[[i], 1].item())
                index = int(table_species_corpora[table_species_corpora["SPECIE"] == specie]["index"])
                # Create association table
                # Prepare data
                table_species_corpora["weights"] = weights[:, index].tolist()
                cooc = table_coocurences[table_coocurences["MESH"] == mesh][["index", "COOC"]]
                data = pd.merge(table_species_corpora, cooc, on = "index", how = "left").fillna(0)
                # Forget data
                data.loc[data["index"] == index, ["TOTAL_PMID_SPECIE", "COOC"]] = [0, 0]
                # Launch analysis
                MeSH_info = table_mesh_corpora[table_mesh_corpora["MESH"] == mesh]
                p = float(MeSH_info["P"])
                r = computation(index, data, p, float(MeSH_info["alpha_prior"]), float(MeSH_info["beta_prior"]), seq = 0.0001, plot = False)
                # fill with results
                validation_set.iloc[i, 2:6] = list(r)
            validation_set.to_csv("data/Validation/neg_classic/validation_" + str(alpha) + "_" + str(sample_size) + ".csv", index = False)

# Anomalie test

if False:
    weaks = pd.read_csv("data/Validation/Anomalie_detection/weakness_set.csv")
    weaks = pd.concat([weaks, pd.DataFrame(columns = ["Mean", "CDF", "Log2FC", "priorCDFratio"])])
    for i in range(0, len(weaks.index)):
        # get row info
        specie = str(weaks.iloc[[i], 0].item())
        mesh = str(weaks.iloc[[i], 1].item())
        index = int(table_species_corpora[table_species_corpora["SPECIE"] == specie]["index"])
        table_species_corpora["weights"] = weights[:, index].tolist()
        cooc = table_coocurences[table_coocurences["MESH"] == mesh][["index", "COOC"]]
        data = pd.merge(table_species_corpora, cooc, on = "index", how = "left").fillna(0)
        MeSH_info = table_mesh_corpora[table_mesh_corpora["MESH"] == mesh]
        p = float(MeSH_info["P"])
        r = computation(index, data, p, float(MeSH_info["alpha_prior"]), float(MeSH_info["beta_prior"]), seq = 0.0001, plot = False)
        weaks.iloc[i, 2:6] = list(r)
    weaks.to_csv("data/Validation/Anomalie_detection/set_" + str(alpha) + "_" + str(sample_size) + ".csv", index = False)

# START TEST
if False:
    index = int(table_species_corpora[table_species_corpora["SPECIE"] == specie]["index"])
    table_species_corpora.insert(2, "weights", weights[:, index].tolist())
    cooc = table_coocurences[table_coocurences["MESH"] == mesh][["index", "COOC"]]
    data = pd.merge(table_species_corpora, cooc, on = "index", how = "left").fillna(0)
    # Forget data
    # data.loc[data["index"] == index, ["TOTAL_PMID_SPECIE", "COOC"]] = [0, 0]
    MeSH_info = table_mesh_corpora[table_mesh_corpora["MESH"] == mesh]
    p = float(MeSH_info["P"])
    print(p)
    print(str(MeSH_info["alpha_prior"]))
    print(str(MeSH_info["beta_prior"]))
    r = computation(index, data, p, float(MeSH_info["alpha_prior"]), float(MeSH_info["beta_prior"]), seq = 0.0001, plot = True)
    print(r)
# END TEST

if False:
    index = int(table_species_corpora[table_species_corpora["SPECIE"] == specie]["index"])
    r2 = specie_mesh(index, table_coocurences, table_species_corpora, weights, table_mesh_corpora)
    r2.to_csv("data/" + specie + "_" + str(alpha) + ".csv", index = False)


# plt.plot(prior_test.x, prior_test.f)
# plt.show()















# print("> Import MeSH corpora sizes ... ", end = '')
# table_mesh_corpora = import_table(args.mesh_corpora_path)
