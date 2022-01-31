#!/usr/bin/python3

import pandas as pd
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description='give arguments to filter script')
parser.add_argument("-l", nargs = 1, required=True, help='log_file')
parser.add_argument('-s', nargs = 1, required=True, help='start time file')
parser.add_argument('-o', nargs = 1, required=True, help='output_dir')
parser.add_argument('-u', nargs = 1, required=True, help='upfront filter')
parser.add_argument('-a', nargs = 1, required=True, help='file with all circRNAs')

args = parser.parse_args()

log_file = pd.read_csv(open(args.l[0]), sep = "\t")
start_time = open(args.s[0])
out_dir = args.o[0]
up_f = args.u[0]
all_circ = open(args.a[0])


summary = open(out_dir + "/summary_run.txt", "a")

#total_nr_circ = log_file.shape[0]
total_nr_circ = 0
all_circ_dict = {}

for line in all_circ:
	total_nr_circ += 1
	all_circ_dict[line.rstrip().split('\t')[0]] = line.rstrip().split('\t')[1:]

primer_found = log_file['primer_found'].sum()
circ_short = total_nr_circ - log_file.shape[0]
primer_no_design = total_nr_circ - log_file['design'].sum() - circ_short
primer_not_found = total_nr_circ - primer_found - primer_no_design - circ_short
total_nr_primers = log_file['total_primer_pairs'].sum()
primers_passed = log_file['passed'].sum()
failed_spec = log_file['failed_spec'].sum()
failed_snp = log_file['failed_SNP'].sum()
failed_str_temp = log_file['failed_sec_str_temp'].sum()
failed_str_amp = log_file['failed_sec_str_amp'].sum()

primer_found_perc = 100 * primer_found/total_nr_circ
circ_short_perc = 100 * circ_short/total_nr_circ
primer_not_found_perc = 100 * primer_not_found/total_nr_circ
primer_no_design_perc = 100 * primer_no_design/total_nr_circ
primers_passed_perc = 100 * primers_passed/ total_nr_primers
failed_spec_perc = 100 * failed_spec / total_nr_primers
failed_snp_perc = 100 * failed_snp / total_nr_primers
failed_str_temp_perc = 100 * failed_str_temp / total_nr_primers
failed_str_amp_perc = 100 * failed_str_amp / total_nr_primers


if up_f == "yes":
	summary.write("{total_nr_circ} circRNA(s) were given as input\n\tfor {primer_found} circRNA(s) ({primer_found_perc:.1f} %), primer pair(s) were designed and filtered succesfully \n\tfor {circ_short} circRNA(s) ({circ_short_perc:.1f} %), the (spliced) template sequence was to short and no primer pair(s) could be designed by primer3\n\tfor {primer_no_design} circRNA(s) ({primer_no_design_perc:.1f} %), no primer pair(s) could be designed by primer3\n\tfor {primer_not_found} circRNA(s) ({primer_not_found_perc:.1f} %), none of the primer pair(s) passed all filters\n\nin total\n\t{total_nr_primers} primer pair(s) were designed and tested\n\t{primers_passed} ({primers_passed_perc:.1f} %) primer pairs passed all filters\n\t{failed_spec} ({failed_spec_perc:.1f} %) primer pair(s) failed the specificity filter\n\t{failed_str_amp} ({failed_str_amp_perc:.1f} %) primer pair(s) failed the secondary structure of the amplicon filter\n\nsee log file for more details per circRNA\n\n".format(total_nr_circ = total_nr_circ, primer_found = primer_found, primer_found_perc = primer_found_perc, primer_no_design = primer_no_design, primer_no_design_perc = primer_no_design_perc, total_nr_primers = total_nr_primers, primers_passed = primers_passed, primers_passed_perc = primers_passed_perc, failed_spec = failed_spec, failed_spec_perc = failed_spec_perc, failed_str_amp = failed_str_amp, failed_str_amp_perc = failed_str_amp_perc, primer_not_found = primer_not_found, primer_not_found_perc = primer_not_found_perc, circ_short = circ_short, circ_short_perc = circ_short_perc))
else:
	summary.write("{total_nr_circ} circRNA(s) were given as input\n\tfor {primer_found} circRNA(s) ({primer_found_perc:.1f} %), primer pair(s) were designed and filtered succesfully \n\tfor {circ_short} circRNA(s) ({circ_short_perc:.1f} %), the (spliced) template sequence was to short and no primer pair(s) could be designed by primer3\n\tfor {primer_no_design} circRNA(s) ({primer_no_design_perc:.1f} %), no primer pair(s) could be designed by primer3\n\tfor {16} circRNA(s) ({17:.1f} %), none of the primer pair(s) passed all filters\n\nin total\n\t{total_nr_primers} primer pair(s) were designed and tested\n\t{primers_passed} ({primers_passed_perc:.1f} %) primer pair(s) passed all filters\n\t{failed_spec} ({failed_spec_perc:.1f} %) primer pair(s) failed the specificity filter\n\t{failed_snp} ({failed_snp_perc:.1f} %) primer pair(s) failed the SNP filter\n\t{failed_str_temp} ({failed_str_temp_perc:.1f} %) primer pair(s) failed the secondary structure of the template filter\n\t{failed_str_amp} ({failed_str_amp_perc:.1f} %) primer pair(s) failed the secondary structure of the amplicon filter\n\nsee log file for more details per circRNA\n\n".format(total_nr_circ = total_nr_circ, primer_found = primer_found, primer_found_perc = primer_found_perc, primer_no_design = primer_no_design, primer_no_design_perc = primer_no_design_perc, total_nr_primers = total_nr_primers, primers_passed = primers_passed, primers_passed_perc = primers_passed_perc, failed_spec = failed_spec, failed_spec_perc = failed_spec_perc, failed_snp = failed_snp, failed_snp_perc = failed_snp_perc, failed_str_temp = failed_str_temp, failed_str_temp_perc = failed_str_temp_perc, failed_str_amp = failed_str_amp, failed_str_amp_perc = failed_str_amp_perc, primer_not_found, primer_not_found_perc, circ_short = circ_short, circ_short_perc = circ_short_perc))

start_str = start_time.read().rstrip()
start = datetime.strptime(start_str, "%d/%m/%Y %H:%M:%S")
end = datetime.now()
end_str = end.strftime("%d/%m/%Y %H:%M:%S")

diff = int((end - start).total_seconds())

summary.write("run started at " + start_str + " (UTC)\nrun ended at " + end_str + " (UTC)\n" + "total run time: " + str(diff) + " seconds\n")

summary.close()


# also add circRNAs that are not in log and filtered primers
log_file = open(args.l[0])

log_file_list = []
for line in log_file:
	log_file_list.append(line.split('\t')[0])

log_file.close()
log_file = open(args.l[0], 'a')
filtered = open('filtered_primers.txt', 'a')

for ID in all_circ_dict:
	if ID not in log_file_list:
		log_file.write(ID + ('\tthis (spliced) circRNA template sequence is too short to design primers for\n'))
		filtered.write(ID + ('\tthis (spliced) circRNA template sequence is too short to design primers for\n'))

log_file.close()
filtered.close()



