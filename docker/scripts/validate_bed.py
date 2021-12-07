#!/usr/bin/python3


import argparse
import sys
import os

parser = argparse.ArgumentParser(description='give arguments to validate bed script')
parser.add_argument('-i', nargs=1, required=True, help='input circRNA bed file, 0-based')
parser.add_argument('-c', nargs=1, required=True, help='chrom file')
args = parser.parse_args()

input_bed = open(args.i[0])
chrom_file = open(args.c[0])

chrom_sizes = {}

try:
	for line in chrom_file:
		c = line.rstrip().split("\t")
		chrom_sizes[c[0]] = int(c[1])
except:
	print("Error parson chrom.sizes file")
	raise

line_number = 0

error = ""

for line in input_bed:
	line_number+=1
	l = line.rstrip().split("\t")

	### check chrom
	chrom = l[0]
	try:
		chromosome_size = chrom_sizes[chrom]
	except:
		error = error + ("Error in bed input file: faulty chromosome name at line %d \n" % (line_number))
		continue

	### check start
	start = l[1]
	try:
		start = int(start)
	except:
		error = error + ("Error in bed input file: parsing start coordinate at line %d \n" % (line_number))
		continue

	if start >= chromosome_size:
		error = error + ("Error in bed input file: start coordinate %d larger than chromosome size %d at line %d \n" % (start, chromosome_size, line_number))
		continue


	### check end
	end = l[2]
	try:
		end = int(end)
	except:
		error = error + ("Error in bed input file: parsing end coordinate at line %d \n" % (line_number))
		continue
	
	if end > chromosome_size:
		error = error + ("Error in bed input file: end coordinate %d larger than chromosome size %d at line %d \n" % (end, chromosome_size, line_number))
		continue

	### check start smaller than start
	if start >= end:
		error = error + ("Error in bed input file: start coordinate %d larger than (or equal to) end coordinate %d at line %d \n" % (start, end, line_number))
		continue

if error != "":
	raise SystemExit(error)
