#!/usr/bin/python3

# # import all libraries
import os
from Bio import Entrez, SeqIO
import argparse
import math

# # get info on BSJ seq
parser = argparse.ArgumentParser(description='give arguments to main primer_xc script')
parser.add_argument('-n', nargs=1, required=True, help='the nr of nucleotides surrounding the BSJ at each side', metavar='length')
parser.add_argument('-i', nargs=1, required=True, help='input circRNA bed file, 0-based')
parser.add_argument('-s', nargs=1, required=True, help='spliced? yes or no')


args = parser.parse_args()
length = int(args.n[0])
input_bed = open(args.i[0])
splice = args.s[0]

# # read first (and only) line of input bed file
circRNA = input_bed.read()

# retrieve chrom start end info
chrom = circRNA.split()[0]
start = int(circRNA.split()[1]) + 1 # change to 1-based system
end = int(circRNA.split()[2])

# # check if circ length is smaller than required BSJ template size

circ_length = end - start

if circ_length < 2 * length:
    length = math.floor(circ_length / 2)

# # retrieve circRNA BSJ sequence

fasta = open("fasta_in.txt", 'w')

#if splice == "yes": 
# ## retrieve left side of BSJ

fasta.write(chrom.replace("chr", "") + ":" + str(end-length) + "-" + str(end) + "\n")

# ## retrieve right side of BSJ

fasta.write(chrom.replace("chr", "") + ":" + str(start) + "-" + str(start + length) + "\n")

# else:

    