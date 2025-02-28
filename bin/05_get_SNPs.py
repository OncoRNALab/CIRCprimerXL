#!/usr/bin/env python3

# # import all libraries
import os
import argparse

parser = argparse.ArgumentParser(description='give arguments to main snp script')
parser.add_argument('-i', nargs=1, required=True, help='input circRNA bed file, 0-based')
parser.add_argument('-u', nargs=1, required=True, help='SNP url', metavar='SNP url')
parser.add_argument('-f', nargs=1, required=True, help='file that keeps track of which sequences are used')
parser.add_argument('-t', nargs=1, required=True, help='file that keeps track of which sequences are used on both ends')

args = parser.parse_args()
input_bed = args.i[0]
SNP_url = args.u[0]

input_file = open(input_bed)
circRNA = input_file.read()

# retreive chrom start end info and add nr of bases bases

chrom = circRNA.split()[0]
start = int(circRNA.split()[1])
end = int(circRNA.split()[2])
circ_ID = circRNA.split()[3]
length = int(circRNA.split()[4])

output_file = open("output_SNP_" + circ_ID + ".bed", "w")


if SNP_url == 'off':
	avoid_range = '[]'

else:

	# retreive SNPs in that region
	os.system("bigBedToBed " + SNP_url + " -chrom="+ chrom + " -start=" + str(start) + " -end=" + str(end) + " " + circ_ID + "_tmp.bed")
	
	
	# get info on which sequences were used
	fasta = open(args.f[0])
	fasta_track = open(args.t[0])

	poss_end = []
	poss_start = []

	for seq, seq_type in zip(fasta, fasta_track):
		
		chrom, start, end = seq.rstrip().replace(":", "-").split("-")
		start = int(start) - 1  	# 1-based to 0-based
		end = int(end)

		if seq_type[0:4] == '>end':
			poss_end = list(range(start,end)) + poss_end

		elif seq_type[0:6] == '>start':
			poss_start = poss_start + list(range(start,end))

	# this variable is a list of all circRNA positions in the right order for template seq
	poss = poss_end + poss_start	
	

	# make list of SNPs in that region
	SNPs = open(circ_ID + "_tmp.bed")
	
	SNP_all = []
	avoid_range = []
	
	for SNP in SNPs:
		SNP_pos = int(SNP.split()[1])
		SNP_length = int(SNP.split()[2]) - int(SNP.split()[1])
	
		for pos in range(SNP_length):
			SNP_all.append(SNP_pos)
			SNP_pos += 1


	# retrieve those SNPs that are actually in the template sequence
	
	for snp in SNP_all:
		if snp in poss:
			avoid_range.append(poss.index(snp))  # index can be used as positions in poss are in correct order (they are 0-based and index as well so no correction needed)


output_file.write(str(avoid_range))

output_file.close()

