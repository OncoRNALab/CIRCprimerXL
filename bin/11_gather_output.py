#!/usr/bin/env python3

import argparse
import os

parser = argparse.ArgumentParser(description='give arguments to main primer_xc script')
parser.add_argument('-i', nargs=1, required=True, help='input circRNA bed file, 0-based')
parser.add_argument('-a', nargs=1, required=True, help='annotation of start and end circ pos')
args = parser.parse_args()
input_file = args.i[0]

all_primers = open(input_file)
circ_ID = str(input_file).split('_')[3]
output = open('selected_primers_' + circ_ID + '.txt', "w")

# get splice info
splice_info = open(args.a[0])
circ_id, info = splice_info.readline().rstrip().split('**')


if os.path.getsize(input_file) == 0:
	circ_ID = str(input_file).split('_')[3]
	chrom = str(input_file).split('_')[4]
	start = str(input_file).split('_')[5]
	end = str(input_file).split('_')[6]

	output.write(circ_ID + '\t' + chrom + '\t' + start + '\t' + end + "\tprimer3 was not able to design primers for this circRNA; try less strict settings\n")

else: 
	primer_found = "no"
	# filter primers according to the filter string
	for primer in all_primers:
		filter_str = primer.split()[-1]
		if (primer_found == "no") and (filter_str == "PASS"):
			output.write(primer.rstrip() + "\t" + info.lstrip('/') + "\n")
			primer_found = 'yes'
	if primer_found == "no":
		output.write(('\t'.join(primer.split()[0:4])) + "\tno primer pair passed all filters for this circRNA\n")


output.close()

		