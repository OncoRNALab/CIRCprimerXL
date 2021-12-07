#!/usr/bin/python3

# # import all libraries
import os
import argparse

parser = argparse.ArgumentParser(description='give arguments to main snp script')
parser.add_argument('-i', nargs=1, required=True, help='input circRNA bed file, 0-based')
parser.add_argument('-u', nargs=1, required=True, help='SNP url', metavar='SNP url')
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
	
	
	start_poss = list(range(start, start + length))
	end_poss = list(range(end - length, end))
	
	
	SNPs = open(circ_ID + "_tmp.bed")
	
	SNP_all = []
	avoid_range = []
	
	for SNP in SNPs:
		SNP_pos = int(SNP.split()[1])
		SNP_length = int(SNP.split()[2]) - int(SNP.split()[1])
	
		for pos in range(SNP_length):
			SNP_all.append(SNP_pos)
			SNP_pos += 1
	
	for snp in SNP_all:
		if snp in start_poss:
			avoid_range.append(length + snp - start)
		if snp in end_poss:
			avoid_range.append(length - (end - snp))


output_file.write(str(avoid_range))

output_file.close()

