#!/usr/bin/python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs=1, required=True, help='input gtf file')
parser.add_argument('-o', nargs=1, required=True, help='output bed file')

args = parser.parse_args()

gtf = open(args.i[0])
bed = open(args.o[0], 'w')

for line in gtf:
	chrom, havana, annotation, start, end, point, strand, point2, info = line.split('\t')

	if annotation == "exon":
		gene_id, gene_version, transcript_id, transcript_version, exon_number, skip = info.split(';', 5)
		
		transcript_id = transcript_id.lstrip().rstrip().replace('transcript_id ', '').replace('"', '')
		transcript_version = transcript_version.lstrip().rstrip().replace('transcript_version ', '').replace('"', '')
		exon_number = exon_number.lstrip().rstrip().replace('exon_number ', '').replace('"', '')

		bed.write("chr" + chrom + '\t' + str(int(start) - 1 ) + '\t' + end + '\t' + transcript_id + '_' + transcript_version + '_exon_' + exon_number + '\n')
