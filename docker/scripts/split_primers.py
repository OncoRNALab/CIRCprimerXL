#!/usr/bin/python3

import argparse

parser = argparse.ArgumentParser(description='give arguments to main primer_xc script')
parser.add_argument('-i', nargs=1, required=True, help='input primer file')
args = parser.parse_args()
input_primers = args.i[0]

primer_in = open(input_primers)

primers = {}

circ_info_keys = ("SEQUENCE_ID", "SEQUENCE_TEMPLATE", "SEQUENCE_TARGET")

# read all info into dictionary
for line in primer_in:
	key, value = line.split("=")
	value = value.rstrip()
	primers[key] = value

template = primers["SEQUENCE_TEMPLATE"]
circRNA = primers["SEQUENCE_ID"]
circ_ID, chrom, start, end = circRNA.split("_")
nr_p_out = primers["PRIMER_LEFT_NUM_RETURNED"]

primer_in.close()


# read general info into file
general_info = open("general_primer_design_info_" + circ_ID + ".txt", "a")
for info in primers:
	if "_NUM_" in info or "_EXPLAIN" in info or any(x in info for x in circ_info_keys):
		general_info.write(info + '=' + str(primers[info]) +'\n')

general_info.close()

# make file for bowtie
primer_file = open("primer_spec_input_" + circ_ID + ".txt", "a")

# make general file with list primers
all_primers = open("all_primers_" + circ_ID + ".txt", 'w')
all_amplicon = open("amplicon_folding_input_" + circ_ID + ".txt", 'w')

for primer_index in range(int(nr_p_out)):


	FWD = primers[("PRIMER_LEFT_" + str(primer_index) + "_SEQUENCE")]
	FWD_qual = len(FWD) * "I"
	REV = primers[("PRIMER_RIGHT_" + str(primer_index) + "_SEQUENCE")]
	REV_qual = len(REV) * "I"

	PRIMER_LEFT_TM = primers[("PRIMER_LEFT_" + str(primer_index) + "_TM")]
	PRIMER_RIGHT_TM = primers[("PRIMER_RIGHT_" + str(primer_index) + "_TM")]
	PRIMER_LEFT_GC_PERCENT = primers[("PRIMER_LEFT_" + str(primer_index) + "_GC_PERCENT")]
	PRIMER_RIGHT_GC_PERCENT = primers[("PRIMER_RIGHT_" + str(primer_index) + "_GC_PERCENT")]

	# bowtie input file
	# write FWD + REV
	primer_file.write(circ_ID + "_primer_" + str(primer_index) + "_FWD_REV" + "\t")
	primer_file.write(FWD + "\t" + FWD_qual + "\t" + REV + "\t" + REV_qual + "\n")

	# write  REV + FWD
	primer_file.write(circ_ID + "_primer_" + str(primer_index) + "_REV_FWD" + "\t")
	primer_file.write(REV + "\t" + REV_qual + "\t" + FWD + "\t" + FWD_qual + "\n")


	# write FWD + FWD
	primer_file.write(circ_ID + "_primer_" + str(primer_index) + "_FWD_FWD" + "\t")
	primer_file.write(FWD + "\t" + FWD_qual + "\t" + FWD + "\t" + FWD_qual + "\n")

	# write REV + REV
	primer_file.write(circ_ID + "_primer_" + str(primer_index) + "_REV_REV" + "\t")
	primer_file.write(REV + "\t" + REV_qual + "\t" + REV + "\t" + REV_qual + "\n")

	# get amplicon and make file for NUPACK
	FWD_pos, FWD_len = primers['PRIMER_LEFT_'+ str(primer_index)].split(",")
	REV_pos, REV_len = primers['PRIMER_RIGHT_'+ str(primer_index)].split(",")
	amplicon = template[int(FWD_pos):int(REV_pos) + 1]


	all_amplicon.write("> amplicon_" + circ_ID + "_primer" + str(primer_index) + "_" + amplicon + "\n")


	# general primer file (for filtering)
	all_primers.write(circ_ID + "\t" + chrom + "\t" + start + "\t" + end + '\t' + str(primer_index) + '\t' + FWD + '\t' + REV + '\t' + 
	FWD_pos + '\t' + FWD_len + '\t' + REV_pos +'\t' + REV_len + '\t' + PRIMER_LEFT_TM + '\t' + PRIMER_RIGHT_TM + '\t' + 
	PRIMER_LEFT_GC_PERCENT + '\t' + PRIMER_RIGHT_GC_PERCENT + '\t' + amplicon + '\n')


primer_file.close()
all_primers.close()
all_amplicon.close()



