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

N = 8877780

g = import_metabolic_network(args.g_path)

print("> Import species corpora sizes ... ", end = '')
table_species_corpora = import_and_map_indexes(args.specie_corpora_path, g)
print("Ok")

print("> Import species-MeSH co-occurences ... ", end = '')
table_coocurences = import_and_map_indexes(args.specie_mesh_path, g)
print("Ok")

print("> Import MeSH corpora sizes ... ", end = '')
table_mesh_corpora = import_table(args.mesh_corpora_path)

# Compute MeSH probabilities
table_mesh_corpora["P"] = table_mesh_corpora["TOTAL_PMID_MESH"]/N

# Compute prior parameters:
mesh_priors = table_mesh_corpora["TOTAL_PMID_MESH"].apply(estimate_prior_distribution_mesh)
mesh_priors = pd.DataFrame(mesh_priors.tolist(), columns = ["alpha_prior", "beta_prior"])
table_mesh_corpora = pd.concat([table_mesh_corpora, mesh_priors], axis = 1)
# table_mesh_corpora = table_mesh_corpora.head(100)
print("Ok")


mesh = "D018312"
specie = "M_tststerones"
# index = 1115
alpha = 0

probabilities = propagation_volume(g, alpha = alpha)

# cc = (100 * probabilities.FOT).round(3)
# cc.to_csv("FOT_" + str(alpha) + ".csv")

if False:
    validation_set = pd.read_csv("data/validation_set_associations.csv")
    # add results columns
    validation_set = pd.concat([validation_set, pd.DataFrame(columns = ["Mean", "CDF", "Log2FC", "priorCDFratio"])])
    # Iter over associations
    for i in range(0, len(validation_set.index)):
        # get row info
        specie = str(validation_set.iloc[[i], 0].item())
        mesh = str(validation_set.iloc[[i], 1].item())
        index = int(table_species_corpora[table_species_corpora["SPECIE"] == specie]["index"])
        # Create association table
        # Prepare data
        table_species_corpora["weights"] = probabilities.FOT.iloc[:, index].tolist()
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
    validation_set.to_csv("data/validation_out.csv", index = False)

# START TEST
if True:
    index = int(table_species_corpora[table_species_corpora["SPECIE"] == specie]["index"])
    table_species_corpora.insert(2, "weights", probabilities.FOT.iloc[:, index].tolist())
    cooc = table_coocurences[table_coocurences["MESH"] == mesh][["index", "COOC"]]
    data = pd.merge(table_species_corpora, cooc, on = "index", how = "left").fillna(0)
    # Forget data
    data.loc[data["index"] == index, ["TOTAL_PMID_SPECIE", "COOC"]] = [0, 0]
    data.to_csv("data/toto.csv", index = False)
    MeSH_info = table_mesh_corpora[table_mesh_corpora["MESH"] == mesh]
    p = float(MeSH_info["P"])
    r = computation(index, data, p, float(MeSH_info["alpha_prior"]), float(MeSH_info["beta_prior"]), seq = 0.0001, plot = True)
    print(r)
# END TEST

if False:
    r2 = specie_mesh(index, table_coocurences, table_species_corpora, probabilities.FOT, table_mesh_corpora)
    r2.to_csv("data/tests/test_full2.csv", index = False)


# plt.plot(prior_test.x, prior_test.f)
# plt.show()

