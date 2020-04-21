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
- `-o, --output OUTFILE`
  - Output file as a GenBank flat file
 
### insect_mito_stats.py
Python script to calculate nucleotide composition stats (AT%, AT-skew, and GC-skew) from an input file in GenBank for FASTA format. Outputs a tab-delimited table.

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
- `-o, --output OUTFILE`
  - Output file as a tab-delimited table

