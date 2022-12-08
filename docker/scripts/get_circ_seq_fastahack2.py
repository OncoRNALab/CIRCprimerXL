#!/usr/bin/python3

# # import all libraries
import os
from Bio import Entrez, SeqIO
import argparse
import math

# # get info on BSJ seq
parser = argparse.ArgumentParser(description='give arguments to main primer_xc script')
parser.add_argument('-i', nargs=1, required=True, help='input circRNA bed file, 0-based')
parser.add_argument('-p', nargs=1, required=True, help='nr of primers')
parser.add_argument('-n', nargs=1, required=True, help='nr of difference between primers')

parser.add_argument('-a', nargs=1, required=True, help='min TM')
parser.add_argument('-b', nargs=1, required=True, help='max TM')
parser.add_argument('-c', nargs=1, required=True, help='opt TM')
parser.add_argument('-d', nargs=1, required=True, help='TM diff')
parser.add_argument('-e', nargs=1, required=True, help='min GC')
parser.add_argument('-f', nargs=1, required=True, help='max GC')
parser.add_argument('-g', nargs=1, required=True, help='opt GC')
parser.add_argument('-j', nargs=1, required=True, help='min amp length')
parser.add_argument('-k', nargs=1, required=True, help='max amp length')



args = parser.parse_args()
input_bed = open(args.i[0])
nr = args.p[0]
diff = args.n[0]

# # read first (and only) line of input bed file
circRNA = input_bed.read()

# retrieve chrom start end info
chrom = circRNA.split()[0]
start = int(circRNA.split()[1]) + 1 # change to 1-based system
end = int(circRNA.split()[2])
circ_ID = circRNA.split()[3]


# retrieve sequence
fasta = open("fasta_out.txt")
fasta_track = open("fasta_track.txt")

right = ""
left = ""

for seq, seq_type in zip(fasta, fasta_track):
	if seq_type[0:4] == '>end':
		left = seq.rstrip() + left
	elif seq_type[0:6] == '>start':
		right = right + seq.rstrip()

# ## paste both sides of BSJ together
sequence = left + right

length = len(sequence) / 2

# ## get max amp length => by user or depending on temp_l
amp_max = args.k[0]
if int(amp_max) == 0:
	amp_max = str(int(length + 30))

# ## make a txt file for primer design (sequence based)

output = open("input_primer3_" + circ_ID + ".txt", "w")
output.write("SEQUENCE_ID=" + circ_ID + "_" + chrom + "_" + str(start-1) + "_" + str(end) + "\n")
output.write("SEQUENCE_TEMPLATE=" + sequence + "\n")
output.write("SEQUENCE_TARGET=" + str(int(length-1)) + ",1\nPRIMER_NUM_RETURN=" + nr + "\n")
output.write("PRIMER_MIN_THREE_PRIME_DISTANCE=" + diff + '\n')
output.write("PRIMER_PRODUCT_SIZE_RANGE=" + args.j[0] + '-' + amp_max + '\n')
output.write("PRIMER_MIN_TM=" + args.a[0] + '\n')
output.write("PRIMER_MAX_TM=" + args.b[0]+ '\n')
output.write("PRIMER_OPT_TM=" + args.c[0]+ '\n')
output.write("PRIMER_PAIR_MAX_DIFF_TM=" + args.d[0]+ '\n')
output.write("PRIMER_MIN_GC="  + args.e[0]+ '\n')
output.write("PRIMER_MAX_GC=" + args.f[0]+ '\n')
output.write("PRIMER_OPT_GC_PERCENT="  + args.g[0]+ '\n')
output.close()


# ## make fasta file folding (sequence based)

output = open("input_NUPACK_" + circ_ID + ".fasta", "w")
output.write("> " + circ_ID + "\n")
output.write(sequence + "\n")
output.close()


# ## make a bed file for SNPs (0-based => use original bed file annotation)

output = open("input_SNP_" + circ_ID + ".bed", 'w')
output.write(circRNA + '\t' + str(int(length)) + '\n')
output.close()

# ## make file for filtering

output = open("input_filter_" + circ_ID + ".bed", 'w')
output.write(circRNA + '\n')
output.close()


