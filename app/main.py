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


mesh = "D000138"
# specie = "M_acorn"
index = 42
alpha = 0

probabilities = propagation_volume(g, alpha = alpha)

# cc = (100 * probabilities.FOT).round(3)
# cc.to_csv("FOT_" + str(alpha) + ".csv")


# START TEST
if False:
    table_species_corpora.insert(2, "weights", probabilities.FOT.iloc[:, index].tolist())
    cooc = table_coocurences[table_coocurences["MESH"] == mesh][["index", "COOC"]]
    data = pd.merge(table_species_corpora, cooc, on = "index", how = "left").fillna(0)
    MeSH_info = table_mesh_corpora[table_mesh_corpora["MESH"] == mesh]
    p = float(MeSH_info["P"])
    r = computation(index, data, p, float(MeSH_info["alpha_prior"]), float(MeSH_info["beta_prior"]), seq = 0.0001)
    print(r)
# END TEST

if True:
    r2 = specie_mesh(42, table_coocurences, table_species_corpora, probabilities.FOT, table_mesh_corpora)
    r2.to_csv("test.csv", index = False)


# plt.plot(prior_test.x, prior_test.f)
# plt.show()

