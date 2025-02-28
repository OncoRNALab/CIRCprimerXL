#!/usr/bin/env python3


import argparse

parser = argparse.ArgumentParser(description='give arguments to filter script')
parser.add_argument("-i", nargs = 1, required=True, help='original primer3 file')
parser.add_argument("-a", nargs = 1, required=True, help='template folding to avoid')
parser.add_argument("-s", nargs=1, required=True, help='snp out')
parser.add_argument("-f", nargs=1, required=True, help='shoud this step be included? no/yes/snp/str')
parser.add_argument('-l', nargs=1, required=True, help='the nr of nucleotides surrounding the BSJ at each side')



args = parser.parse_args()

length = int(args.l[0])
primer3_file = open(args.i[0])

circ_ID, chrom, start, end = primer3_file.readline().replace("SEQUENCE_ID=", "").split("_")
start = int(start)
end = int(end)


primer3_file.close()

filter_on = args.f[0]
primer3_file = open(args.i[0])

primer3_file_new = open("primer3_file_" + circ_ID + ".txt", 'w')


# copy complete file
for line in primer3_file:
	primer3_file_new.write(line)
	# save length info
	if line[0:17] == "SEQUENCE_TEMPLATE":
		length = len(line.split("=")[1]) / 2

avoid_range = []

# if this filter is on
if filter_on == 'yes' or filter_on == "str":
	folding_file = open(args.a[0])

	# get folding info
	skip = folding_file.readline()
	skip = folding_file.readline()
	skip = folding_file.readline()
	fold_temp_avoid = folding_file.readline()
	folding_file.close()
	

	# put folding info into format for primer 3 (space-separated: position,length)
	# only if there are positions to avoid
	if fold_temp_avoid != '[]':
		fold_temp_avoid = fold_temp_avoid.replace('[', "").replace(']', '').split(', ')
		fold_temp_avoid = [int(x) for x in fold_temp_avoid]

		begin_region = fold_temp_avoid[0]
		index = 0
		length = 1

		for position in fold_temp_avoid:
			if index + 1 < len(fold_temp_avoid):
				if position == fold_temp_avoid[index + 1] - 1:
					length += 1
					index += 1
				else:
					avoid_range.append([begin_region, length])
					begin_region = fold_temp_avoid[index + 1]
					index += 1
					length = 1
			else:
				avoid_range.append([begin_region, length])
				length = 1

if filter_on == 'yes' or filter_on == "snp":

	SNPs = open(args.s[0]).readline()

	if SNPs != '[]':
		SNPs = SNPs.replace('[', "").replace(']', '').split(', ')
		SNPs = [ int(x) for x in SNPs ]

		for snp in SNPs:
			avoid_range.append([snp, 1])

if filter_on != "no":
		avoid_range_str = "SEQUENCE_EXCLUDED_REGION="

		for element in avoid_range:
			avoid_range_str = avoid_range_str + str(element[0]) + "," + str(element[1]) + ' '

		primer3_file_new.write(avoid_range_str.rstrip() + '\n')


# add last line with '=' again
primer3_file_new.write("=\n")

primer3_file.close()
primer3_file_new.close()


