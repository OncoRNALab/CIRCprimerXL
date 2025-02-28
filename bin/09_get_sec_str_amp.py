#!/usr/bin/env python3

from nupack import *
import argparse

parser = argparse.ArgumentParser(description='give arguments to main primer_xc script')
parser.add_argument('-i', nargs=1, required=True, help='nupack input file')
args = parser.parse_args()
input_file = args.i[0]

sequence_file = open(input_file)
circ_ID = str(input_file).split('_')[3]
circ_ID = circ_ID.replace(".txt", "")
output = open('output_NUPACK_amp_' + circ_ID + '.txt', 'a')

for line in sequence_file:
	tmp, tmp2, primer_ID, seq = line.split("_")
	sequence = []
	sequence.append(seq.rstrip())

	# make model
	my_model = Model(material='DNA', celsius=60, sodium=0.05, magnesium=0.003)

	# MFE proxy structure and energy
	my_mfe = mfe(strands=sequence, model=my_model)
	delta_g = float(my_mfe[0].energy)
	structure = str(my_mfe[0].structure)

	# make string from regions to avoid
	str_not_ok = []
	index = 0

	for sign in structure:
		if sign == '(' or sign == ')':
			str_not_ok.append(index)
		index += 1

	# Print out components of the result for the given complex
	
	output.write(circ_ID + '_' + primer_ID +'\t' + str(delta_g) + '\t' + str(str_not_ok) + '\t' + structure + '\n')
	
output.close()
	