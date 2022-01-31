#!/usr/bin/python3

# # import all libraries
import os
from Bio import Entrez, SeqIO
import argparse
import math
import pybedtools
import re

# # get info on BSJ seq
parser = argparse.ArgumentParser(description='give arguments to main primer_xc script')
parser.add_argument('-n', nargs=1, required=True, help='the nr of nucleotides surrounding the BSJ at each side', metavar='length')
parser.add_argument('-i', nargs=1, required=True, help='input circRNA bed file, 0-based')
parser.add_argument('-s', nargs=1, required=True, help='spliced? yes or no')
parser.add_argument('-e', nargs=1, required=True, help='bed file with canonical exons')
parser.add_argument('-t', nargs=1, required=True, help='list with canonical transcripts')



args = parser.parse_args()
temp_length = int(args.n[0])
input_bed = open(args.i[0])
splice = args.s[0]

# # read first (and only) line of input bed file
circRNA = input_bed.read()

# retrieve chrom start end info
circ_chrom = circRNA.split()[0]
circ_start = int(circRNA.split()[1]) #+ 1 # change to 1-based system
circ_end = int(circRNA.split()[2])
circ_ID = circRNA.split()[3]


# add parameters for circ annotation
annotation = open('annotation_' + circ_ID + '.txt', 'w')
annotation.write(circ_ID + "**")

# # retrieve circRNA BSJ sequence

fasta = open("fasta_in.txt", 'w')
fasta_track = open("fasta_track.txt", 'w')

# to keep track of nr of exons used (and order!)
fasta_start = 1 
fasta_end = 1

########## STEP 1
# ## unspliced
if splice == "no":

    # check if circ length is smaller than required BSJ template size
    circ_length = circ_end - circ_start
    
    if circ_length < 2 * temp_length:
        temp_length = math.floor(circ_length / 2)

    # retrieve left side of BSJ
    fasta_track.write(">end_" + str(fasta_end) + '\n')
    fasta.write(circ_chrom.replace("chr", "") + ":" + str(circ_end - temp_length + 1) + "-" + str(circ_end) + "\n") # +1 to change to 1-based system
    
    # retrieve right side of BSJ
    fasta_track.write(">start_" + str(fasta_start) + '\n')
    fasta.write(circ_chrom.replace("chr", "") + ":" + str(circ_start + 1) + "-" + str(circ_start + temp_length) + "\n") # +1 to change to 1-based system
    annotation.write('not_computed' + '\t')
    annotation.write('not_computed' + '\t')
    annotation.write('unspliced')


########## STEP 1
# ## spliced
else:
    
    # check if circ length is smaller than required BSJ template size
    circ_length = circ_end - circ_start
    
    if circ_length < 2 * temp_length:
        temp_length = math.floor(circ_length / 2)

    # set up bed file with canonical exon positions
    exons = args.e[0]

    ## see if END pos BSJ is in exon
    # make file with theoretical left side of BSJ
    end_bed = open("end.bed", "w")
    end_bed.write(circ_chrom + "\t" + str(circ_end - temp_length) + "\t" + str(circ_end) + "\n") # 1-based system always needs -1 to when substracting positions 
    end_bed.close()
    
    # check if END falls into an exon
    # make both files bedtools format
    exons_bed = pybedtools.BedTool(exons)
    end_bed = pybedtools.BedTool('end.bed')
    # intersect both files
    match_end = exons_bed.intersect(end_bed)
    # write results to new file
    match_file_end = open('match_file_end.bed', 'w')
    match_file_end.write(str(match_end))
    match_file_end.close()
    
    # make dict of results
    match_file_end = open('match_file_end.bed')
    match_dict_end_size = {}
    match_dict_end_pos = {}
    match_ex_tr_end = {}
    
    for line in match_file_end:
        match_c, match_s, match_e, trans = line.split('\t')
        match_dict_end_size[trans] = int(match_e) - int(match_s) + 1
        match_dict_end_pos[trans] = match_c + ':' + match_s + '-' + match_e
        match_ex_tr_end[trans.split("_")[0]] = trans

    ## see if START pos BSJ is in exon
    # make file with theoretical right side of BSJ
    start_bed = open("start.bed", "w")
    start_bed.write(circ_chrom + "\t" + str(circ_start) + "\t" + str(circ_start + temp_length) + "\n")
    start_bed.close()
    
    # check if START falls into an exon
    # make both files bedtools format
    exons_bed = pybedtools.BedTool(exons)
    start_bed = pybedtools.BedTool('start.bed')
    # intersect both files
    match_start = exons_bed.intersect(start_bed)
    # write results to new file
    match_file_start = open('match_file_start.bed', 'w')
    match_file_start.write(str(match_start))
    match_file_start.close()
    
    # make dict of results
    match_file_start = open('match_file_start.bed')
    match_dict_start_size = {}
    match_dict_start_pos = {}
    match_ex_tr_start = {}
    
    for line in match_file_start:
        match_c, match_s, match_e, trans = line.split('\t')
        match_dict_start_size[trans] = int(match_e) - int(match_s) + 1
        match_dict_start_pos[trans] = match_c + ':' + match_s + '-' + match_e
        match_ex_tr_start[trans.split("_")[0]] = trans


    ########## STEP 2
    ################# if one of them not in exon => do unspliced
    if len(match_dict_end_size) == 0 or len(match_dict_start_size) == 0:
        fasta_track.write(">end_" + str(fasta_end) + '\n')
        fasta.write(circ_chrom.replace("chr", "") + ":" + str(circ_end - temp_length + 1) + "-" + str(circ_end) + "\n") # change to 1-based system
        fasta_track.write(">start_" + str(fasta_start) + '\n')
        fasta.write(circ_chrom.replace("chr", "") + ":" + str(circ_start + 1) + "-" + str(circ_start + temp_length) + "\n") # change to 1-based system
        # add info to annotation file
        if len(match_dict_end_size) == 0:
            annotation.write('intronic' + '\t')
        else:
            annotation.write('exonic' + '\t')
        if len(match_dict_start_size) == 0:
            annotation.write('intronic' + '\t')
        else:
            annotation.write('exonic' + '\t')

        annotation.write('unspliced')
    
    
    ########## STEP 2
    else:
        
        # make list of ENST that are common between start and end
        trans_overlap = list(set(match_ex_tr_start.keys()).intersection(set(match_ex_tr_end.keys())))

        ########## STEP 3
        # if there are none => only use the two current largest exons (might be long enough)
        if len(trans_overlap) == 0:

            ### END of BSJ
            # get the largest exon 
            max_exon_end_length = 0
            
            for match, match_pos in match_dict_end_size.items():
                if match_pos > max_exon_end_length:
                    max_exon_end = match
                    max_exon_end_length = match_pos
            
            # this variable describes how many nts are currently in the template length on the end side
            temp_length_end = circ_end - int(match_dict_end_pos[max_exon_end].replace(':', '-').split('-')[1])

            if temp_length_end >= temp_length:
                # if it's long enough, the template sequence can be generated as usual
                fasta_track.write(">end_" + str(fasta_end) + '\n')
                fasta.write(circ_chrom.replace("chr", "") + ":" + str(circ_end - temp_length + 1) + "-" + str(circ_end) + "\n") # change to 1-based system
                annotation.write('exon_' + match_dict_end_pos[max_exon_end] + '\t')

            # if not, it needs the whole exon can be used
            else:
                # print that exon pos to fasta file
                fasta_track.write(">end_" + str(fasta_end) + '\n')
                fasta.write(circ_chrom.replace("chr", "") + ":" + str(circ_end - temp_length_end + 1) + "-" + str(circ_end) + "\n") # change to 1-based system
                

            # START of BSJ
            # get the largest exon 
            max_exon_start_length = 0
            
            for match, match_pos in match_dict_start_size.items():
                if match_pos > max_exon_start_length:
                    max_exon_start = match
                    max_exon_start_length = match_pos
            
            # this variable describes how many nts are currently in the template length on the start side
            temp_length_start = int(match_dict_start_pos[max_exon_start].replace(':', '-').split('-')[2]) - circ_start
        
            if temp_length_start >= temp_length:
                # if it's long enough, the template sequence can be generated as usual
                fasta_track.write(">start_" + str(fasta_start) + '\n')
                fasta.write(circ_chrom.replace("chr", "") + ":" + str(circ_start + 1) + "-" + str(circ_start + temp_length) + "\n") # change to 1-based system
                annotation.write('exon_' + match_dict_start_pos[max_exon_start] + '\t')

            # if not, the whole exon can be used
            else:

                # print that exon pos to fasta file
                fasta_track.write(">start_" + str(fasta_start) + '\n')
                fasta.write(circ_chrom.replace("chr", "") + ":" + str(circ_start + 1) + "-" + str(circ_start + temp_length_start) + "\n") # change to 1-based system


        ########## STEP 3
        # use common transcript (trans_overlap)
        else:

            # read in file with canonical transcripts if a file is provided
            trans_list = []
            selected_trans = ""
            if args.t[0] != "":
                trans_file = open(args.t[0])
                for line in trans_file:
                    trans_list.append(line.rstrip())

                # check if common transcript is canonical
                if len(list(set(trans_overlap).intersection(set(trans_list)))) > 0:
                    selected_trans = list(set(trans_overlap).intersection(set(trans_list)))[0]
        
            # if no file provided or no canonical trans in list => select other transcript
            ########## STEP 4
            if len(selected_trans) == 0:
                selected_trans = trans_overlap[0]
            

            # get info about the current exons
            ### END of BSJ
            new_exon_end = match_ex_tr_end[selected_trans]

            # this variable describes how many nts are currently in the template length on the end side
            temp_length_end = circ_end - int(match_dict_end_pos[new_exon_end].replace(':', '-').split('-')[1])
            # this variable describes the left most position of the start template side (to avoid overlap with start template)
            new_end = circ_end - temp_length_end
         
            ### START of BSJ
            new_exon_start = match_ex_tr_start[selected_trans]

            # this variable describes how many nts are currently in the template length on the start side
            temp_length_start = int(match_dict_start_pos[new_exon_start].replace(':', '-').split('-')[2]) - circ_start
            # this variable describes the right most position of the start template side (to avoid overlap with end template)
            new_start = circ_start + temp_length_start

    ############################## add END sequence
            ########## STEP 5
            if temp_length_end >= temp_length:
                # if it's long enough, the template sequence can be generated as usual
                fasta_track.write(">end_" + str(fasta_end) + '\n')
                fasta.write(circ_chrom.replace("chr", "") + ":" + str(circ_end - temp_length + 1) + "-" + str(circ_end) + "\n") # change to 1-based system
                annotation.write('exon_' + match_dict_end_pos[new_exon_end] + '\t')

            # if not, it needs to be recalculated
            ########## STEP 5
            else:

                # print that exon pos to fasta file
                fasta_track.write(">end_" + str(fasta_end) + '\n')
                fasta.write(circ_chrom.replace("chr", "") + ":" + str(circ_end - temp_length_end + 1) + "-" + str(circ_end) + "\n") # change to 1-based system
                new_end = int(circ_end - temp_length_end)
                annotation_end = 'exon_' + match_dict_end_pos[new_exon_end]

                extra_exon_end = 0  # var needed to augment exon nr if more than one extra exon is required


                # while loop to add as many exons as needed
                while temp_length_end < temp_length:

                    # look for previous exon and save their end pos in separate dict
                    # split depending on organism (c elengans has differen annotation
                    if new_exon_end.count('_') == 2:
                        trans_nr, skip, exon_nr = new_exon_end.split("_", 2)

                    else:
                        trans_nr, trans_v, skip, exon_nr = new_exon_end.split("_", 3)

                    exon_nr = int(exon_nr)
                    previous_exon = trans_nr + '_exon_' + str(exon_nr - 1 - extra_exon_end)

                    # open file with canonical exons
                    exons_file = open(exons)

                    # retrieve previous exon
                    new_exons_end_pos = {}

                    for line in exons_file:
                        if re.search(previous_exon, line):
                            new_c, new_s, new_e, new_trans = line.split("\t")
                            # add new exons to the dict with most right pos
                            new_exons_end_pos[new_c.replace("chr", "") + ":" + new_s + "-" + new_e] = int(new_e)

                    exons_file.close()

                    ####### STEP 6
                    # if no new exon is found
                    if len(new_exons_end_pos) == 0:
                        temp_length_end = 10000

                    else:

                        # select the exon that's most close by based on pos dict
                        selected_exon_end = max(new_exons_end_pos)
                        selected_c, selected_s, selected_e = selected_exon_end.replace(':', '-').split('-')


                        ########## STEP 7
                        # make sure that exon is not the one that's already used for the START pos (this is the case for 2-exon circs)
                        # if start of new selected exon is bigger than end of exon at START pos, then it can be used            
                        if int(selected_exon_end.replace(':', '-').split('-')[1]) > new_start:
                            # get correct pos of that exon
                            # if new exon is long enough, take the part you need
                            if int(selected_e) - int(selected_s) >= temp_length - temp_length_end:
                                selected_seq = selected_c + ":" + str(
                                    int(selected_e) - (temp_length - temp_length_end)) + "-" + selected_e
                                selected_seq_1 = selected_c + ":" + str(
                                    int(selected_e) - (temp_length - temp_length_end) + 1) + "-" + selected_e # change to 1-based system

                                new_end = int(selected_e) - (temp_length - temp_length_end)
                                temp_length_end = temp_length


                            # if not, take the whole exon
                            else:
                                selected_seq = selected_c + ":" + selected_s + "-" + selected_e
                                selected_seq_1 = selected_c + ":" + str(int(selected_s) + 1) + "-" + selected_e # change to 1-based system
                                temp_length_end = temp_length_end + (int(selected_e) - int(selected_s))
                                new_end = int(selected_s)
                                extra_exon_end += 1

                           # write it to fasta file
                            fasta_end += 1
                            fasta_track.write(">end_" + str(fasta_end) + '\n')
                            fasta.write(selected_seq_1 + "\n") # change to 1-based system
                            annotation_end = annotation_end + '/exon_' + selected_seq


                        ########## STEP 7
                        else:
                            # set temp length end to 10.000 so while loop stops
                            temp_length_end = 10000
                
                annotation.write(annotation_end + '\t')


    ################################# add START sequence
            
            ########## STEP 5
            if temp_length_start >= temp_length:
                # if it's long enough, the template sequence can be generated as usual
                fasta_track.write(">start_" + str(fasta_start) + '\n')
                fasta.write(circ_chrom.replace("chr", "") + ":" + str(circ_start + 1) + "-" + str(circ_start + temp_length) + "\n") # change to 1-based system
                annotation.write('exon_' + match_dict_start_pos[new_exon_start] + '\t')
            
            ########## STEP 5
            # if not, it needs to be recalculated
            else:

                fasta_track.write(">start_" + str(fasta_start) + '\n')
                fasta.write(circ_chrom.replace("chr", "") + ":" + str(circ_start + 1) + "-" + str(circ_start + temp_length_start) + "\n") # change to 1-based system
                new_start = int(circ_start + temp_length_start)
                annotation_start = 'exon_' + match_dict_start_pos[new_exon_start]

                extra_exon_start = 0 # var needed to augment exon nr if more than one extra exon is required

                # while loop to add as many exons as needed
                while temp_length_start < temp_length:

                    # look for next exon and save their start pos in separate dict
                    # split depending on organism (c elengans has differen annotation
                    if new_exon_end.count('_') == 2:
                        trans_nr, skip, exon_nr = new_exon_start.split("_", 2)

                    else:
                        trans_nr, trans_v, skip, exon_nr = new_exon_start.split("_", 3)
                    exon_nr = int(exon_nr)
                    next_exon = trans_nr + '_exon_' + str(exon_nr + 1 + extra_exon_start)

                    # open file with canonical exons
                    exons_file = open(exons)

                    # retrieve next exon
                    new_exons_start_pos = {}

                    for line in exons_file:
                        if re.search(next_exon, line):
                            new_c, new_s, new_e, new_trans = line.split("\t")
                            # add new exons to the dict with most left pos
                            new_exons_start_pos[new_c.replace("chr", "") + ":" + new_s + "-" + new_e] = int(new_s)

                    exons_file.close()

                    ####### STEP 6
                    # if no new exon is found
                    if len(new_exons_start_pos) == 0:
                        temp_length_start = 10000

                    else:
                        # select the exon that's most close by based on pos dict
                        selected_exon_start = min(new_exons_start_pos)
                        selected_c, selected_s, selected_e = selected_exon_start.replace(':', '-').split('-')
      
                        ########## STEP 7                    
                        # make sure that exon is not the one that's already used for the END pos (this is the case for 2-exon circs)
                        # if end of new selected exon is smaller than start of exon at END pos, then it can be used

                        if int(selected_exon_start.replace(':', '-').split('-')[2]) < new_end:
                            # get correct pos of that exon
                            # if new exon is long enough, take the part you need
                            if int(selected_e) - int(selected_s) >= temp_length - temp_length_start:
                                selected_seq = selected_c + ":" + selected_s + "-" + str(
                                    int(selected_s) + (temp_length - temp_length_start))
                                selected_seq_1 = selected_c + ":" + str(int(selected_s) + 1) + "-" + str(
                                    int(selected_s) + (temp_length - temp_length_start)) # change to 1-based system

                                new_start = int(selected_s) + (temp_length - temp_length_start)
                                temp_length_start = temp_length


                            # if not, take the whole exon
                            else:
                                selected_seq = selected_c + ":" + selected_s + "-" + selected_e
                                selected_seq_1 = selected_c + ":" + str(int(selected_s) + 1) + "-" + selected_e # change to 1-based system
                                temp_length_start = temp_length_start + (int(selected_e) - int(selected_s))
                                new_start = int(selected_e)
                                extra_exon_start += 1

                            # write it to fasta file
                            fasta_start += 1
                            fasta_track.write(">start_" + str(fasta_start) + '\n')
                            fasta.write(selected_seq_1 + "\n") # change to 1-based system
                            annotation_start =  annotation_start + '/exon_' + selected_seq
                        
                        ########## STEP 7
                        # set temp length start to 10.000 so while loop stops
                        else:
                            temp_length_start = 10000

                annotation.write(annotation_start + '\t')

            annotation.write('spliced')

annotation.write('\n')

fasta.close()
fasta_track.close()
annotation.close()
