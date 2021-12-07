#!/usr/bin/python3

import argparse
import os

parser = argparse.ArgumentParser(description='give arguments to main primer_xc script')
parser.add_argument('-i', nargs=1, required=True, help='input circRNA bed file, 0-based')
args = parser.parse_args()
input_file = args.i[0]

all_primers = open(input_file)
circ_ID = str(input_file).split('_')[2]
output = open('selected_primers_' + circ_ID + '.txt', "w")

# output.write("circ_ID\tchr\tstart\tend\tprimer_ID\tFWD_primer\tREV_primer\tFWD_pos\tFWD_length\tREV_pos\tREV_length\tFWD_TM\tREV_TM\tFWD_GC\tREV_GC\tamplicon\n")

if os.path.getsize(input_file) == 0:
	circ_ID = str(input_file).split('_')[3]
	chrom = str(input_file).split('_')[4]
	start = str(input_file).split('_')[5]
	end = str(input_file).split('_')[6]

	output.write(circ_ID + '\t' + chrom + '\t' + start + '\t' + end + "\tprimer3 was not able to design primers for this circRNA; try less strict settings\n")

else: 
	primer_found = "no"
	for primer in all_primers:
		filter_str = primer.split()[-1]
		if (primer_found == "no") and (filter_str == "PASS"):
			output.write(primer)
			primer_found = 'yes'
	if primer_found == "no":
		output.write(('\t'.join(primer.split()[0:4])) + "\tno primer pair passed all filters for this circRNA\n")


output.close()

		