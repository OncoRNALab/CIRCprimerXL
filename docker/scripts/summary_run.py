#!/usr/bin/python3

import pandas as pd
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description='give arguments to filter script')
parser.add_argument("-l", nargs = 1, required=True, help='log_file')
parser.add_argument('-s', nargs = 1, required=True, help='start time file')
parser.add_argument('-o', nargs = 1, required=True, help='output_dir')


args = parser.parse_args()

log_file = pd.read_csv(open(args.l[0]), sep = "\t")
start_time = open(args.s[0])
out_dir = args.o[0]

summary = open(out_dir + "/summary_run.txt", "a")

total_nr_circ = log_file.shape[0]
primer_found = log_file['primer_found'].sum()
primer_no_design = total_nr_circ - log_file['design'].sum() 
primer_not_found = total_nr_circ - primer_found - primer_no_design
total_nr_primers = log_file['total_primer_pairs'].sum()
primers_passed = log_file['passed'].sum()
failed_spec = log_file['failed_spec'].sum()
failed_snp = log_file['failed_SNP'].sum()
failed_str_temp = log_file['failed_sec_str_temp'].sum()
failed_str_amp = log_file['failed_sec_str_amp'].sum()

primer_found_perc = 100 * primer_found/total_nr_circ
primer_not_found_perc = 100 * primer_not_found/total_nr_circ
primer_no_design_perc = 100 * primer_no_design/total_nr_circ
primers_passed_perc = 100 * primers_passed/ total_nr_primers
failed_spec_perc = 100 * failed_spec / total_nr_primers
failed_snp_perc = 100 * failed_snp / total_nr_primers
failed_str_temp_perc = 100 * failed_str_temp / total_nr_primers
failed_str_amp_perc = 100 * failed_str_amp / total_nr_primers

summary.write("{0} circRNAs were given as input\n\tfor {1} circRNAs ({2:.1f} %), primer pairs were designed and filtered succesfully \n\tfor {3} circRNAs ({4:.1f} %), no primer pairs could be designed by primer3\n\tfor {16} circRNAs ({17:.1f} %), none of the primer pairs passed all filters\n\nin total\n\t{5} primer pairs were designed and tested\n\t{6} ({7:.1f} %) primer pairs passed all filters\n\t{8} ({9:.1f} %) primer pairs failed the specicifity filter\n\t{10} ({11:.1f} %) primer pairs failed the SNP filter\n\t{12} ({13:.1f} %) primer pairs failed the secundary structure of the template filter\n\t{14} ({15:.1f} %) primer pairs failed the secundary structure of the amplicon filter\n\nsee log file for more details per circRNA\n\n".format(total_nr_circ, primer_found, primer_found_perc, primer_no_design, primer_no_design_perc, total_nr_primers, primers_passed, primers_passed_perc, failed_spec, failed_spec_perc, failed_snp, failed_snp_perc, failed_str_temp, failed_str_temp_perc, failed_str_amp, failed_str_amp_perc, primer_not_found, primer_not_found_perc))

start_str = start_time.read().rstrip()
start = datetime.strptime(start_str, "%d/%m/%Y %H:%M:%S")
end = datetime.now()
end_str = end.strftime("%d/%m/%Y %H:%M:%S")

diff = int((end - start).total_seconds())

summary.write("run started at " + start_str + " (UTC)\nrun ended at " + end_str + " (UTC)\n" + "total run time: " + str(diff) + " seconds\n")

summary.close()