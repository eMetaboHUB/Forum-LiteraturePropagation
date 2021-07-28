# Propagation

## Environment

You need conda installed

```bash
conda env create --file env/propagation.yml
```
Then, use :

```bash
conda activate propagation
```



On prend en compte 2 types de graphs :
    - Undirected: les probas seront calculées à partir de la matrice d'adjascence
    - Directed: les noeuds doivent contenir un attribut "Weight" qui correspond aux probabilité de transition entre les noeuds

Pour les valeurs retourner, on renvoie également TOTAL_PMID_SPECIE et COOC mais pas TOTAL_PMID_MESH, car étudiant les mentions ce n'est pas la même valeurs que dans l'analyse classique.