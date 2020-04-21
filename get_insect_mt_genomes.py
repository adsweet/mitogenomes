#!/usr/bin/env python3

# Written by Andrew D. Sweet
# 10 April 2020
# Usage: python3 get_insect_mt_genomes.py -e <email address> -o <output file>

from Bio import SeqIO
from Bio import Entrez
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--email","-e", help="Email address for Entrez", required=True)
parser.add_argument("--output","-o", help="Output file", required=True)
args = parser.parse_args()

Entrez.email = args.email
search = "(txid6960[Organism:exp] AND ((minicircle[TI] OR minichromosome[TI]) AND mitochondrion[TI)) OR (txid6960[Organism:exp] AND mitochondrion, complete genome[TI] AND 10000:20000[SLEN])"
handle = Entrez.esearch(db="nucleotide", term=search, rettype="gb", retmode="xml")
record = Entrez.read(handle)

count = int(record["Count"])
print(count)
handle2 = Entrez.esearch(db="nucleotide", term=search, rettype="gb", retmode="xml", retmax=count)
record2 = Entrez.read(handle2)

id_list = record2["IdList"]
post = Entrez.epost("nucleotide", id=",".join(id_list))
search_results = Entrez.read(post)

webenv = search_results["WebEnv"]
key = search_results["QueryKey"]

out = open(args.output, "w")
gbf = Entrez.efetch(db="nucleotide", webenv=webenv, query_key=key, rettype="gb", retmode="text", retmax=count)
gbf_record = gbf.read()
gbf.close()
out.write(gbf_record)
out.close()

handle.close()
handle2.close()









