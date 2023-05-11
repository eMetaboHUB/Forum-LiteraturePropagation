# Suggesting disease associations for overlooked metabolites using literature from metabolic neighbours

<div style="text-align: justify">
In human health research, metabolic signatures extracted from metabolomics data are a strong-added value for stratifying patients and identifying biomarkers. Nevertheless, one of the main challenges is to interpret and relate these lists of discriminant metabolites to pathological mechanisms. This task requires experts to combine their knowledge with information extracted from databases and the scientific literature.


However, we show that a large fraction of metabolites are rarely or never mentioned in the literature. Consequently, these overlooked metabolites are often set aside and the interpretation of metabolic signatures is restricted to a subset of the significant metabolites. To suggest potential pathological phenotypes related to these understudied metabolites, we extend the 'guilt by association' principle to literature information by using a Bayesian framework. With this approach, we suggest more than 35,000 associations between 1,047 overlooked metabolites and 3,288 diseases (or disease families). All these newly inferred associations are freely available on the FORUM ftp server (See information below).
</div>

## Availability

All computed relations are available on the FORUM ftp server.

- host: sftp://ftp.semantic-metabolomics.org
- user: forum
- password: Forum2021Cov!

They can be accessed using any FTP client (like FileZilla) or the command line interface sftp. For downloading all the suggested associations, use for instance:

```bash
sftp forum@ftp.semantic-metabolomics.org:/Propagation/Human1_1.7/2021/global.csv
```

The directory contains: 

- *global.csv*: all the computed associations.
- *predictions.csv*: the selected associations with: LogOdds > 2, Log2FC > 1, Entropy > 1.
- *specie_names*: species labels in the Human1 1.7 metabolic network.
- *mesh_labels*: MeSH descriptor labels.
- *overlooked_metabolites_list.tsv*: the list of the 2113 metabolites considered as overlooked in the analysis carried out in the article. This file is also stored in the data directory of this repository.

In the RDF directory, you can find the RDF graphs related to these associations.

```bash
sftp forum@ftp.semantic-metabolomics.org:/Propagation/Human1_1.7/2021/RDF/*.ttl.gz
```

Supplementary data and html reports for predictions:

The **data.tar.gz** archive contains all the contributors tables for each of the 35585 extracted and provided associations. Reciprocally, the **report.tar.gz** contains all the associated html reports. All the reports are also available online through static html files. See for example for the relationship between [5alpha-androstane-3,17-dione and Polycystic Ovary Syndrome](https://forum-static-files.semantic-metabolomics.fr/1.1/report/report_M_m01064c_D011085.html). See more info [here](https://forum-static-files.semantic-metabolomics.fr/).


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

  --undirected: Is the graph undirected ? If the graph is directed, the probability are extracted from the 'Weight' attribute of the edges. In case of a multi-graph, the multiple edges are merged by summing the weights. If the graph in undirected, transition probabilities are computed from the adjacency matrix

  --specie.corpora: Path to the SPECIE corpus file

  --specie.cooc:    Path to the SPECIE-MeSH co-occurences file

  --mesh:   The studied MeSH descriptor identifier (eg. D010661). This argument is incompatible with the 'file' argument. This will return all associations between this MeSH and all species in the metabolic network. If the 'specie' option is also set, only this specific association will be computed.

  --specie: The studied metabolic specie identifier in the network according to the *id_att* value (eg. M_m01778c). This argument is incompatible with the 'file' argument. The program will return all associations between this specie and all MeSHs. If the 'mesh' option is also set, only this association specific association will be computed.

  --file:   Path to a file containing pairs of SPECIE and MESH associations to be computed (format csv: SPECIE, MESH). This argument is incompatible with the 'mesh' and 'specie' arguments.

  --forget: Only the prior from neighbourhood will be used in the computation, observations are hidden.

  --alpha:  The damping factor (alpha). It could be a single value, or several values to test different parametrizations. All provided alpha values will be tested against all provided sample size values. Default = 0.4

  --sample_size:    The sample size parameter. It could be a single value, or several values to test different parametrizations. All provided sample size values will be tested against all provided alpha values. Default = 1000

  --q:  The tolerance threshold for PPR probabilities. If q is negative, the default value is used. The Default = 1/(N - 1)

  --id_att: The name of the vertex attribute containing the SPECIE identifier (eg. M_m0001c) in the graph. These identifiers must match with those provided in the 'SPECIE' column of specie.corpora and specie.cooc files. Default is the 'label' attribute

  --species_name_path:  Path to the file containing species names in a .csv format. The first column should be named 'SPECIE' and contains the species' identifiers (eg. M_m0001c) and the second column should be named 'SPECIE_NAME' and contains the species' chemical names

  --meshs_name_path:    Path to the file containing MeSH names in a .csv format. The first column should be named 'MESH' and contains the MeSH identifiers (eg. D000017) and the second column should be named 'MESH_NAME' and contains the MeSH labels

  --out:    Path to the output directory
```

### Data format details

* ```--specie.corpora```: a two-column table

Eg.

| SPECIE    | TOTAL_PMID_SPECIE |
|-----------|-------------------|
| M_m01773c | 9                 |
| M_m01778c | 215               |
| M_m01784c | 6193              |
| ...       | ...               |

* ```--specie.cooc```: a three-column table

Eg.

| SPECIE    | MESH       | COOC |
|-----------|------------|------|
| M_m01778c | D000079426 | 2    |
| M_m01778c | D000230    | 3    |
| M_m01778c | D000236    | 1    |
| ...       | ...        |      |

* ```--file```: an association file with two column

Eg.

| SPECIE    | MESH    |
|-----------|---------|
| M_m03120c | D005271 |
| M_m03099c | D006973 |
| M_m02756c | D020969 |
| ...       | ...     |

* ```--species_name_path```: a two-column file

Eg.

| SPECIE     | SPECIE_NAME                        |
|------------|------------------------------------|
| M_CE5840_c | 2,3-Epoxy-Alpha-Tocopheryl Quinone |
| M_m02731c  | phosphatidate-LD-PS pool           |
| M_m02864c  | S-(PGJ2)-glutathione               |
| ...        | ...                                |

* ```--meshs_name_path```: a two-column file

Eg.

| MESH       | MESH_NAME              |
|------------|------------------------|
| D000017    | ABO Blood-Group System |
| D000068556 | Interferon beta-1a     |
| D000069463 | Olive Oil              |
| ...        | ...                    |

Data are availbale in the *data* directory

## Output

Predictions are reported in a table format with columns:

### Literature and identifiers

- SPECIE: The metabolic specie identifier

- MESH: The MeSH identifier

- TOTAL_PMID_SPECIE: The total number of articles discussing the metabolic specie (from *specie.corpora*)

- COOC: The number of co-mentions between the specie and the MeSH (from *specie.cooc*)

- Mean: The computed averaged probability that one article mentioning the specie would mention the disease

### Estimators

- CDF: The probability that an article mentioning the specie, would mention the disease less frequently than expected (the marginal probability of mentioning the disease)

- LogOdds: The LogOdds computed from the CDF: (1 - CDF) / (CDF)

- Log2FC: Fold-Change computed from the averaged probability that one article mentioning the compound would mention the disease (MEAN), compared to the expected probability (the marginal probability of mentioning the disease).

**Warnings**: 
  - When a specie **has no** available literature (COOC = 0): CDF, LogOdds, and Log2FC are computed from the **prior** distribution (the neighbouring literaure)

  - When a specie **has** annotated literature (COOC > 0): CDF, LogOdds, and Log2FC are computed from the **posterior** distribution (using neighbouring literature and observations)

  - When a specie has **no literature in its neighbourhood** (NeighborhoodInformation = FALSE), a prior can't be built from the neighbouring literature. In this case, we use the *general* prior (see publication for details) and the predictions are computed from its posterior distribution.

  - When a specie has both no annotated literature and no literature in the neighbourdhood (COOC = 0 & NeighborhoodInformation = FALSE), predictions are directly derived from the *general* prior, but should be discarded.

### Diagnostic values

- priorLogOdds: This indicators is computed only for species that have annotated literature (COOC > 0). It corresponds to the LogOdds computed from the prior distribution (the neighbouring literature), without considering the literature of the targeted specie.

- priorLog2FC: Same as priorLogOdds but for Log2FC

- NeighborhoodInformation: Indicates if there is available literature in the neighbourhood.

- Entropy: Entropy of the contributions in the prior distribution

- CtbAvgDistance: Averaged distance of the contributors

- CtbAvgCorporaSize:	Averaged number of articles for the contributors 

- NbCtb: Number of contributors

- SPECIE_NAME: label of the specie (from species_name_path)

- MESH_NAME: label of the MeSH (from meshs_name_path)

## Examples

### For a specific pair of specie and MeSH (```--specie``` and ```--mesh``` options)

```bash
python app/main.py --graph="data/Human1/1.7/Human-GEM_CarbonSkeletonGraph_noComp_no1C_cpds.gml" --specie.corpora="data/Human1/1.7/species_pmids_Human1_1.7.csv" --specie.cooc="data/Human1/1.7/species_mesh_pmids_Human1_1.7.csv" --specie="M_m01064c"  --mesh="D011085" --out="path/to/out/dir"  --alpha 0.4 --sample_size 1000 --species_name_path="data/Human1/1.7/species_names.csv" --meshs_name_path="data/Human1/1.7/mesh_labels.csv"
```

#### 3 files are outputed: 
  - *specie*\_*mesh*\_*alpha*\_*sample_size*.csv: association results

**Eg:** 

| SPECIE    | MESH    | TOTAL_PMID_SPECIE | COOC | Mean | CDF  | LogOdds | Log2FC | priorLogOdds | priorLog2FC | NeighborhoodInformation | Entropy | CtbAvgDistance | CtbAvgCorporaSize | NbCtb | SPECIE_NAME                  | MESH_NAME                 |
|-----------|---------|-------------------|------|------|------|---------|--------|--------------|-------------|-------------------------|---------|----------------|-------------------|-------|------------------------------|---------------------------|
| M_m01064c | D011085 |             82 | 1 | 0.02 | 0.00 |    6.23 |   3.14 |         5.47 |        3.98 | TRUE                    |    2.29 |           1.58 |          35869.47 | 24 | 5alpha-androstane-3,17-dione | Polycystic Ovary Syndrome |

  - data_*specie*\_*mesh*\_*alpha*\_*sample_size*.csv: details of the contributors

**Eg:** 


| SPECIE    | TOTAL_PMID_SPECIE | PriorWeights     | COOC    | PostWeights | PriorLogOdds | PostLogOdds | PriorLog2FC       | PostLog2FC        | SPECIE_NAME          |
|-----------|-------------------|-------------|---------|--------------------|--------------|-------------|-------------------|-------------------|----------------------|
| M_m01338c |              2348 |     0.29    |  45  |        0.51        |     48.59    |    49.44    |        2.60       |        2.59       | androsterone         |
| M_m02969c |             79421 |     0.29    | 2521 |        0.28        |      inf     |     inf     |        3.75       |        3.75       | testosterone         |
| M_m02968c |             69365 |     0.09    | 1939 |        0.10        |      inf     |     inf     |        3.56       |        3.56       | testosterone sulfate |
| M_m01787c |             93909 |     0.02    | 1102 |        0.04        |      inf     |     inf     |        2.32       |        2.32       | estradiol-17beta     |
| M_m01787c | 93909             | 0.023936701 | 37      | 0.024870212        | -inf         | -inf        | -3.65098133939937 | -3.65222726884335 | estradiol-17beta     |

For each contributor: 

- SPECIE: contributor identifier
- TOTAL_PMID_SPECIE: number of annotated articles
- COOC number of co-mentions with the MeSH
- PriorWeights: weights in the prior distribution
- PostWeights: weights in the posterior distribution
- PriorLogOdds: LogOdds (specific to the contributor) computed only from its literature in the **prior** distribution.
- PostLogOdds: LogOdds (specific to the contributor) computed only from its literature in the **posterior** distribution.
- PriorLog2FC: Log2FC (specific to the contributor) computed only from its literature in the **prior** distribution.
- PostLog2FC: Log2FC (specific to the contributor) computed only from its literature in the **posterior** distribution.
- SPECIE_NAME: name of the contributor

**Warnings**: When a specie has **no** literature, the predictions are only based on the prior distribution and therefore PostLogOdds and PostLog2FC are not returned for the contributors. Same for PostWeights which is not returned if the specie has no literature.

  - report_*specie*\_*mesh*\_*alpha*\_*sample_size*.csv: an html report with figure of distributions, a table and contributor profiles.

* * *

### For only a MeSH descriptor (```--mesh``` option):
  - *mesh*.csv: associations with all the metabolic species in the network

```bash
python app/main.py --graph="data/Human1/1.7/Human-GEM_CarbonSkeletonGraph_noComp_no1C_cpds.gml" --specie.corpora="data/Human1/1.7/species_pmids_Human1_1.7.csv" --specie.cooc="data/Human1/1.7/species_mesh_pmids_Human1_1.7.csv" --specie="M_m01064c" --out="path/to/out/dir"  --alpha 0.4 --sample_size 1000 --species_name_path="data/Human1/1.7/species_names.csv" --meshs_name_path="data/Human1/1.7/mesh_labels.csv"
```

| SPECIE    | TOTAL_PMID_SPECIE | COOC | Mean | CDF  | LogOdds | Log2FC | priorLogOdds | priorLog2FC | NeighborhoodInformation | Entropy | CtbAvgDistance | CtbAvgCorporaSize | NbCtb | SPECIE_NAME                               |
|-----------|-------------------|------|------|------|---------|--------|--------------|-------------|-------------------------|---------|----------------|-------------------|-------|-------------------------------------------|
| M_m02731c |              0 | 0 | 0.00 | 1.00 |   -7.57 |  -2.24 |              |             | VRAI                    |    0.89 |           2.24 |          52940.92 | 11 | phosphatidate-LD-PS pool                  |
| M_m02864c |              0 | 0 | 0.00 | 0.99 |   -4.36 |  -2.32 |              |             | VRAI                    |    2.71 |           1.43 |          66854.47 | 73 | S-(PGJ2)-glutathione                      |
| M_m01666c |            282 | 0 | 0.00 | 0.99 |   -5.17 |  -3.39 |        -4.29 |       -3.13 | VRAI                    |    1.29 |           1.17 |          24152.59 | 30 | deoxyadenosine                            |
| M_m01799c |             46 | 8 | 0.01 | 0.00 |   21.51 |   2.60 |         6.12 |        2.20 | VRAI                    |    0.03 |           1.00 |           1162.25 | 14 | etiocholan-3alpha-ol-17-one 3-glucuronide |


* * *


### For only a metabolic specie (```--specie``` option):
  - *specie*.csv: associations with all the MeSH descriptors provided in the *specie.cooc* file.

```bash
python app/main.py --graph="data/Human1/1.7/Human-GEM_CarbonSkeletonGraph_noComp_no1C_cpds.gml" --specie.corpora="data/Human1/1.7/species_pmids_Human1_1.7.csv" --specie.cooc="data/Human1/1.7/species_mesh_pmids_Human1_1.7.csv" --specie="M_m01064c" --out="path/to/out/dir"  --alpha 0.4 --sample_size 1000 --species_name_path="data/Human1/1.7/species_names.csv" --meshs_name_path="data/Human1/1.7/mesh_labels.csv"
```

| MESH    | TOTAL_PMID_SPECIE | COOC | Mean | CDF  | LogOdds | Log2FC | priorLogOdds | priorLog2FC | NeighborhoodInformation | Entropy | CtbAvgDistance | CtbAvgCorporaSize | NbCtb | MESH_NAME               |
|---------|-------------------|------|------|------|---------|--------|--------------|-------------|-------------------------|---------|----------------|-------------------|-------|-------------------------|
| D011469 |             82 | 9 | 0.06 | 0.00 |   23.47 |   3.44 |         4.16 |        2.75 | VRAI                    |    2.29 |           1.58 |          35869.47 | 24 | Prostatic Diseases      |
| D005834 |             82 | 8 | 0.05 | 0.00 |   21.19 |   3.38 |         5.79 |        2.77 | VRAI                    |    2.29 |           1.58 |          35869.47 | 24 | Genital Neoplasms, Male |
| D011471 |             82 | 8 | 0.05 | 0.00 |   19.82 |   3.43 |         4.17 |        2.61 | VRAI                    |    2.29 |           1.58 |          35869.47 | 24 | Prostatic Neoplasms     |
| D005832 |             82 | 9 | 0.12 | 0.00 |   19.15 |   3.50 |         5.68 |        2.89 | VRAI                    |    2.29 |           1.58 |          35869.47 | 24 | Genital Diseases, Male  |

* * *

### For a file of specific association (```--file``` option)
  - *file_name*\_*alpha*\_*sample_size*.csv: associations results


| SPECIE    | MESH    | TOTAL_PMID_SPECIE | COOC | Mean | CDF  | LogOdds | Log2FC | priorCDF | priorLog2FC | NeighborhoodInformation | Entropy | CtbAvgDistance | CtbAvgCorporaSize | NbCtb |
|-----------|---------|-------------------|------|------|------|---------|--------|----------|-------------|-------------------------|---------|----------------|-------------------|-------|
| M_m03120c | D005271 |                 0 |    0 | 0.00 | 0.69 |   -0.79 |   0.15 |     0.69 |        0.15 | VRAI                    |    2.29 |           1.16 |          14201.39 |    41 |
| M_m03099c | D006973 |                 0 |    0 | 0.01 | 0.96 |   -3.20 |  -0.16 |     0.96 |       -0.16 | VRAI                    |    1.34 |           1.06 |          92486.93 |    35 |
| M_m02756c | D020969 |                 0 |    0 | 0.03 | 0.03 |    3.33 |   0.37 |     0.03 |        0.37 | VRAI                    |    0.29 |           1.03 |           7310.96 |    12 |
| M_m01830c | D007592 |                 0 |    0 | 0.01 | 0.28 |    0.93 |   0.34 |     0.28 |        0.34 | VRAI                    |    0.75 |           1.17 |            247.84 |     5 |


For more details, see the publication.
