# Mitogenomes
Scripts for working with sequence data from mitochondrial genomes

## Scripts
### get_insect_mt_genomes.py
Python script to pull insect mitochondrial genomes from GenBank using several esearch terms.

```
Usage: python3 get_insect_mt_genomes.py -e <email adress> -o <output file>
```
#### Arguments
- `-h, --help`
  - Lists arguments
- `-e, --email`
  - An email address for pulling GenBank data using Entrez
- `-o, --output`
  - Output file as a GenBank flat file
 
### insect_mito_stats.py
Python script to calculate nucleotide composition stats (AT%, AT-skew, and GC-skew) from an input file in GenBank or FASTA format. Outputs a tab-delimited table.

```
Usage: python3 insect_mito_stats.py -i <input> -f <input file format> -o <output>
```
#### Arguments
- `-h, --help`
  - Lists arguments
- `-i, --input`
  - Input file
- `-f, --format`
  - File format of the input file; "genbank" or "fasta"
- `-o, --output`
  - Output file as a tab-delimited table
### columbicola_mt_interspecies_comp.R
R script to compare mitogenome gene boundaries in different individuals of a single species. The script is set up to compare individuals in species of <i>Columbicola</i> lice, but can be modifid to work with other taxa. 
### louse_at_comp.R
R script to compare nucleotide composition among the mitogenomes of different species of parasitic lice.
### louse_tree_traits.R
R script to reconstruct the evolution and test for phylogenetic signal of mitogenome structure on a phylogenetic tree of parasitic lice.  

