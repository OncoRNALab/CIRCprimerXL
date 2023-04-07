#!/usr/bin/python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs=1, required=True, help='input MANE file')
parser.add_argument('-o', nargs=1, required=True, help='output ENST list file')


mane = open(args.i[0])
out = open(args.i[0], 'w')

out_list = []

for line in mane:
	if line[0] != '#':
		chrom, BestRefSeq, tr_type, start, end, point, strand, point2, info = line.split('\t')

		if tr_type == "exon":
			gene_id, transcript_id, skip = info.split(';', 2)
			out_list.append(transcript_id.lstrip().rstrip().replace('transcript_id', '').replace('"', ''))
		


out_list = set(out_list)

for element in out_list:
	out.write(element[1:] + '\n')

mane.close()
out.close()
