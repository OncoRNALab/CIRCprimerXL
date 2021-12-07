#!/usr/bin/python3


import argparse

parser = argparse.ArgumentParser(description='give arguments to filter script')
parser.add_argument("-i", nargs = 1, required=True, help='original primer3 file')
parser.add_argument("-a", nargs = 1, required=True, help='template folding to avoid')
parser.add_argument("-f", nargs=1, required=True, help='shoud this step be included? yes')

args = parser.parse_args()

filter_on = args.f[0]
primer3_file = open(args.i[0])
folding_file = open(args.a[0])

circ_ID = str(args.i[0]).split("_")[2].replace(".str", "")
primer3_file_new = open("primer3_file_" + circ_ID + ".txt", 'w')


# copy complete file except last line whith '='
for line in primer3_file:
	if len(line.rstrip()) > 1:
		primer3_file_new.write(line)

# if this filter is on
if filter_on == 'yes':

	# get folding info
	skip = folding_file.readline()
	skip = folding_file.readline()
	skip = folding_file.readline()
	fold_temp_avoid = folding_file.readline()
	folding_file.close()
	avoid_range = []

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

	avoid_range_str = "SEQUENCE_EXCLUDED_REGION="

	for element in avoid_range:
		avoid_range_str = avoid_range_str + str(element[0]) + "," + str(element[1]) + ' '

	primer3_file_new.write(avoid_range_str.rstrip() + '\n')


# add last line with '=' again
primer3_file_new.write("=\n")

primer3_file.close()
primer3_file_new.close()


