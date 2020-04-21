#!/usr/bin/env python3

# Written by Andrew D. Sweet
# 10 April 2020
# Usage: python3 insect_mito_stats.py -i <input> -f <input file format> -o <output>


from Bio import SeqIO
from Bio.Seq import Seq
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input","-i", help="Input file format (genbank or fasta)", required=True)
parser.add_argument("--format","-f", help="Input file format", required=True)
parser.add_argument("--output","-o", help="Output file", required=True)
args = parser.parse_args()

#gb_file = "insect_mito_sequences.gb"

with open(args.output, "w") as output:
	if args.format == "genbank":
		output.write("Species" + '\t' + "Order" + '\t' + "Accession" + '\t' + "AT_percent" + '\t' + "GC_percent" + '\t' + "AT_skew" + '\t' + "GC_skew" + '\n')
		for gb_record in SeqIO.parse(open(args.input, "r"), args.format):
			seq = gb_record.seq
			name = gb_record.annotations["organism"]
			acc = gb_record.id
			taxonomy = gb_record.annotations["taxonomy"]
			orders = ["Archaeognatha","Zygentoma","Ephemeroptera","Odonata","Orthoptera","Neuroptera","Phasmatodea","Embioptera","Grylloblattodea","Mantophasmatodea","Plecoptera","Dermaptera","Zoraptera","Mantodea","Blattodea","Psocoptera","Phthiraptera","Thysanoptera","Hemiptera","Hymenoptera","Strepsiptera","Coleoptera","Megaloptera","Raphidioptera","Trichoptera","Lepidoptera","Diptera","Siphonaptera","Mecoptera"]
			order = next((x for x in orders if x in taxonomy), False)
			length = len(gb_record.seq)
			A = seq.count("A")
			T = seq.count("T")
			G = seq.count("G")
			C = seq.count("C")
			AT = (A + T)/length
			GC = (G + C)/length
			AT_skew = (A-T)/(A+T)
			GC_skew = (G-C)/(G+C)
			output.write(name + '\t' + str(order) + '\t' + acc + "\t" + str(AT) + "\t" + str(GC) + "\t" + str(AT_skew) + "\t" + str(GC_skew) + '\n')
	if args.format == "fasta":
		output.write("Species" + '\t' + "AT_percent" + '\t' + "GC_percent" + '\t' + "AT_skew" + '\t' + "GC_skew" + '\n')
		for gb_record in SeqIO.parse(open(args.input, "r"), args.format):
			seq = gb_record.seq
			name = gb_record.id
			length = len(gb_record.seq)
			A = seq.count("A")
			T = seq.count("T")
			G = seq.count("G")
			C = seq.count("C")
			AT = (A + T)/length
			GC = (G + C)/length
			AT_skew = (A-T)/(A+T)
			GC_skew = (G-C)/(G+C)
			output.write(name + '\t' + str(AT) + "\t" + str(GC) + "\t" + str(AT_skew) + "\t" + str(GC_skew) + '\n')
