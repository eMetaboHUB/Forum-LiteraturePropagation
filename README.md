# Propagation

Suggesting relations between metabolic species and disease-related MeSH descriptors by propagating the neighbouring literature.
[ABSTRACT]

## Environment
Requirement: conda

Install env:
```bash
conda env create --file env/propagation.yml
```

Activate env
```bash
conda activate propagation
```

## Usage:

```
python app/main.py --graph="path/to/graph" --specie.corpora="path/to specie-table" --specie.cooc="path/to/specie-mesh-table" --specie="specie"  --mesh="mesh" --out="path/to/output/dir"  --alpha a --sample_size n --species_name_path="/path/to/specie-name-table" --meshs_name_path="path/to/mesh-name-table"


  -h, --help    show this help message and exit

  --graph:  Path to the metabolic network compound graph in .gml format

  --undirected: Is the graph undirected ? If the graph is directed the probability are extracted from the 'Weight' attribute of the edges. In case of a multi-graph, the multiple edges are merged by summing the weights. If the graph in undirected transition probabilities are computed form the adjacency matrix

  --specie.corpora: Path to the SPECIE corpus file

  --specie.cooc:    Path to the SPECIE-MeSH co-occurences file

  --mesh:   The studied MeSH descriptor identifier (eg. D010661). This argument is incompatible with the 'file' argument. This will return all associations between this MeSH and all species in the metabolic network. If the 'specie' option is also set, only this specific association will be computed.

  --specie: The studied specie identifier in the network according to the *id_att* value (eg. M_m01778c). This argument is incompatible with the 'file' argument. The program will return all associations between this specie and all MeSHs. If the 'mesh' option is also set, only this association specific association will be computed.

  --file:   Path to a file containing pairs of SPECIE and MESH associations to be computed (format csv: SPECIE, MESH). This argument is incompatible with the 'mesh' and 'specie' arguments.

  --forget: Only the prior from neighborhood will be used in the computation, observations of treated species are hidden.

  --alpha:  The damping factor (alpha). It could be a single value, or several values to test different parametrizations. All provided alpha values will be tested against all provided sample size values. Default = 0.1

  --sample_size:    The sample size parameter. It could be a single value, or several values to test different parametrizations. All provided sample size values will be tested against all provided alpha values. Default = 100

  --q:  The tolerance threshold for PPR probabilities. If q is negative the default value is used. The Default = 1/(N - 1)

  --id_att: The name of the vertex attribute containing the SPECIE identifier (eg. M_m0001c) in the graph. These identifiers must match with those provided in the 'SPECIE' column of specie.corpora and specie.cooc files. Default is the 'label' attribute

  --species_name_path:  Path to the file containing species names in a .csv format. First column should be named 'SPECIE' and contains the SPECIE identifiers (eg. M_m0001c) and the second column should be named 'SPECIE_NAME' and contains the species' chemical names

  --meshs_name_path:    Path to the file containing species names in a .csv format. First column should be named 'MESH' and contains the MESH identifiers (eg. D000017) and the second column should be named 'MESH_NAME' and contains the MeSH labels

  --out:    Path to the output directory
```

### Data format details

* ```--specie.corpora```: a two column table
Eg.

| SPECIE    | TOTAL_PMID_SPECIE |
|-----------|-------------------|
| M_m01773c | 9                 |
| M_m01778c | 215               |
| M_m01784c | 6193              |
| ...       | ...               |

* ```--specie.cooc```: a three column table
Eg.

| SPECIE    | MESH       | COOC |
|-----------|------------|------|
| M_m01778c | D000079426 | 2    |
| M_m01778c | D000230    | 3    |
| M_m01778c | D000236    | 1    |
| ...       | ...        |      |

* ```--file```: a two column requested association file
Eg.

| SPECIE    | MESH    |
|-----------|---------|
| M_m03120c | D005271 |
| M_m03099c | D006973 |
| M_m02756c | D020969 |
| ...       | ...     |

* ```--species_name_path```: a two column file
Eg.

| SPECIE     | SPECIE_NAME                        |
|------------|------------------------------------|
| M_CE5840_c | 2,3-Epoxy-Alpha-Tocopheryl Quinone |
| M_m02731c  | phosphatidate-LD-PS pool           |
| M_m02864c  | S-(PGJ2)-glutathione               |
| ...        | ...                                |

* ```--meshs_name_path```: a two column file
Eg.

| MESH       | MESH_NAME              |
|------------|------------------------|
| D000017    | ABO Blood-Group System |
| D000068556 | Interferon beta-1a     |
| D000069463 | Olive Oil              |
| ...        | ...                    |

Data are availbale in the *data* directory

## Output

* For a specific pair of specie and MeSH (```--specie``` and ```--mesh``` options), 3 files are outputed: 
   - *specie*\_*mesh*\_*alpha*\_*sample_size*.csv: association results
   - data_*specie*\_*mesh*\_*alpha*\_*sample_size*.csv: details of the prior contributors
  - report_*specie*\_*mesh*\_*alpha*\_*sample_size*.csv: an html report with figure of distributions and contributor profiles

* For only a MeSH descriptor (```--mesh``` option):
  - *mesh*.csv: associations with all the metabolic species in the network

* For only a metabolic specie (```--specie``` option):
  - *specie*.csv: associations with all the MeSH descriptors provided in the *specie.cooc* file.

* For a file of specific association (```--file``` option)
  * *file_name*\_*alpha*\_*sample_size*.csv: associations results



## Example of execution

```bash
python app/main.py --graph="data/Human1/1.7/Human-GEM_CarbonSkeletonGraph_noComp_no1C_cpds.gml" --specie.corpora="data/Human1/1.7/species_pmids_Human1_1.7.csv" --specie.cooc="data/Human1/1.7/species_mesh_pmids_Human1_1.7.csv" --specie="M_34dhpe_c"  --mesh="D010300" --out="notes/results/"  --alpha 0.4 --sample_size 1000 --species_name_path="data/Human1/1.7/species_names.csv" --meshs_name_path="data/Human1/1.7/mesh_labels.csv"
```