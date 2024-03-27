#!/usr/bin/env nextflow

// set default parameters

params.index_bowtie = "$baseDir/assets/GRCh38/index_bowtie"
params.index_bowtie_name = "GRCh38_dna"
params.index_fasta = "$baseDir/assets/GRCh38/index_fastahack"
params.index_fasta_name = "GRCh38.dna.primary_assembly.fa"
params.primer_settings = "$baseDir/assets/primer3plus_settings.txt"
params.chrom_file = "$baseDir/assets/GRCh38/chrom_sizes_GRCh38.txt"
params.input_bed = "example/path"
params.known_exons = "$baseDir/assets/GRCh38/known_exons_GRCh38.bed"
params.list_ENST = "$baseDir/assets/GRCh38/ENST_list_GRCh38.txt"

params.splice = 'yes'
params.primer3_diff = 1
params.primer3_nr = 20
params.min_tm = 58
params.max_tm = 60
params.opt_tm = 59
params.diff_tm = 2
params.min_gc = 30
params.max_gc = 80
params.opt_gc = 50
params.amp_min = 50
params.amp_max = 0 // this param is set to 0, so that it can be adjusted depending on temp_l if the user does not supply amp_max
params.temp_l = 150
params.snp_filter = 'strict'
params.spec_filter = 'strict'
params.snp_url = 'http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153Common.bb'
params.upfront_filter = "yes"

params.help = false

// required parameters
input_bed = file(params.input_bed)
chrom_file = file(params.chrom_file)
index_bowtie = file(params.index_bowtie)
index_fasta = file(params.index_fasta)
known_exons = file(params.known_exons)
list_ENST = file(params.list_ENST)

// help message

def helpMessage() {
	log.info"""
	
	Usage:
	
	The typical command for running the pipeline is as follows:
	nextflow run CIRCprimerXL.nf -profile singularity
	
	Mandatory nextflow arguments:
	-profile 		set to 'local' when running locally, set to 'singularity' when running on the HPC

	Mandatory pipeline arguments:
	--input_bed		path to input file with circRNAs in bed format (0-based annotation)
	--index_bowtie		path to bowtie genome index directory
	--index_bowtie_name	the basename of the Bowtie index to be searched (the name of any of the index files up to but not including the final .1.ebwt / .rev.1.ebwt / ...)
	--index_fasta		path to directory that contains the fastahack genome and index file
	--index_fasta_name	the name of the fastahack genome file


	Optional pipeline arguments:
	--splice			when set to 'yes' the input sequence will be spliced, when set to 'no' the input sequence will be unspliced
	--primer_settings	path to file with primer3plus settings (see primer3 manual)
	--chrom_file		file containing all chromosome sizes (to validate bed file)
	--known_exons		bed file containing exon annotation
	--list_ENST			file containing ENST numbers of canonical transcripts or transcripts of interest (this file can also be left empty)
	--primer3_diff		the minimum number of base pairs between the 3' ends of any two left primers (see also primer3 PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE)
	--primer3_nr		the number of primers designed by primer3; caution: setting this parameter to a large value will increase running time
	--min_tm		minimum melt temperature of the primers (default: 58)
	--max_tm		maximum melt temperature of the primers(default: 60)
	--opt_tm		optimal melt temperature of the primers(default: 59)
	--diff_tm		maximum difference in melt temperature between the primers(default: 2)
	--min_gc		minimum GC contect of the primers (default: 30)
	--max_gc		maximum GC contect of the primers(default: 80)
	--opt_gc		optimal GC contect of the primers(default: 50)
	--amp_min		minimum amplicon length (default: 60)
	--amp_max		maximum amplicon length (default: 0)
	--temp_l		the number of nucleotides on each side of the circRNA BSJ that will be used for the template (example 150 => template of 300 nts in total)
	--spec_filter		when set to 'strict', only 2MM + 3MM are allowed; when set to 'loose', 2MM + 2MM and 2MM + 1MM are also allowed
	--snp_filter		when set to 'strict', no common SNPs are allowed in primer sequence; when set to 'loose', common SNPs are allowed in 5' first half of primer; when set to 'off', no filter is applied
	--snp_url		when using a differente species than human, the correct SNP database url should be provided; alternatively, this paramater can be set to 'off' if no SNP database is available
	--upfront_filter	when set to 'yes', SNPs and secundary structures are avoided before primer design; when set to 'str', secundary structures are avoided before primer design; when set to 'snp', snp are avoided before primer design; when set to 'no', no filtering before primer design is performed
	--output_dir		path to directory where the output files will be saved


	"""
}

if (params.help) {
	helpMessage()
	exit 0
}

// check if all parameters are provided correctly

if (!file(params.input_bed).exists()) {
	exit 1, "Input bed file not found: ${params.input_bed}"}

if (!file(params.output_dir).exists()) {
	exit 1, "Output directory not found: ${params.output_dir}"}

if (!file(index_bowtie).exists()) {
	exit 1, "Index file not found: ${index_bowtie}"}

if (!file(params.primer_settings).exists()) {
	exit 1, "Primer3 settings file not found: ${params.primer_settings}"}

if (!file(params.chrom_file).exists()) {
	exit 1, "Chromosome size file not found: ${params.chrom_file}"}

if (params.splice == "yes"){
	if (!file(params.known_exons).exists()) {
		exit 1, "Known exons file not found: ${params.known_exons}"}}

if (params.splice == "yes" && !(params.list_ENST == "none")) {
	if (!file(params.list_ENST).exists()) {
		exit 1, "ENST list file not found: ${params.list_ENST}"}}

if (!params.primer3_diff.toString().isNumber()){
	exit 1, "Invalid primer3_diff PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE: ${params.primer3_diff}. Valid options: any integer > 0."}

if (params.primer3_diff.toInteger() < 0){
	exit 1, "Invalid primer3_diff PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE: ${params.primer3_diff}. Valid options: any integer > 0."}

if (!params.primer3_nr.toString().isNumber()){
	exit 1, "Invalid primer3_nr PRIMER_NUM_RETURN: ${params.primer3_nr}. Valid options: any integer > 0. Caution: setting this parameter to a large value will increase running time."}
if (params.primer3_nr.toInteger() < 0){
	exit 1, "Invalid primer3_nr PRIMER_NUM_RETURN: ${params.primer3_nr}. Valid options: any integer > 0. Caution: setting this parameter to a large value will increase running time."}

if (!params.temp_l.toString().isNumber()){
	exit 1, "Invalid template length: ${params.temp_l}. Valid options: any integer between 50 and 500."}
if (params.temp_l.toInteger() < 50 || params.temp_l.toInteger() > 500){
	exit 1, "Invalid template length: ${params.temp_l}. Valid options: any integer between 50 and 500."}

if (params.snp_filter != "strict" && params.snp_filter != 'loose' && params.snp_filter != 'off'){
	exit 1, "Invalid SNP filter: ${params.snp_filter}. Valid options: 'strict','loose'."}

if (params.spec_filter != "strict" && params.spec_filter != 'loose'){
	exit 1, "Invalid specificity filter: ${params.spec_filter}. Valid options: 'strict','loose'."}

if (params.splice != "yes" && params.splice != 'no'){
	exit 1, "Invalid splicing option: ${params.splice}. Valid options: 'yes','no'."}


if (params.upfront_filter != "yes" && params.upfront_filter != 'str' && params.upfront_filter != 'snp' && params.upfront_filter != 'no'){
	exit 1, "Invalid SNP filter: ${params.upfront_filter}. Valid options: 'yes','str','snp','no'."}

if (!params.min_tm.toString().isNumber()) {exit 1, " min_tm: ${params.min_tm}. Valid options: any integer > 0."}
if (!params.max_tm.toString().isNumber()) {exit 1, " max_tm: ${params.max_tm}. Valid options: any integer > 0."}
if (!params.opt_tm.toString().isNumber()) {exit 1, " opt_tm: ${params.opt_tm}. Valid options: any integer > 0."}
if (!params.diff_tm.toString().isNumber()) {exit 1, " diff_tm: ${params.diff_tm}. Valid options: any integer > 0."}
if (!params.min_gc.toString().isNumber()) {exit 1, " min_gc: ${params.min_gc}. Valid options: any integer > 0."}
if (!params.max_gc.toString().isNumber()) {exit 1, " max_gc: ${params.max_gc}. Valid options: any integer > 0."}
if (!params.opt_gc.toString().isNumber()) {exit 1, " opt_gc: ${params.opt_gc}. Valid options: any integer > 0."}
if (!params.amp_min.toString().isNumber()) {exit 1, " amp_min: ${params.amp_min}. Valid options: any integer > 0."}
if (!params.amp_max.toString().isNumber()) {exit 1, " amp_max: ${params.amp_max}. Valid options: any integer > 0."}

if (params.min_tm.toInteger() > params.max_tm.toInteger() ) {exit 1, " min_tm and max_tm: max_tm (${params.min_tm}) should be > min_tm (${params.max_tm})"}
if (params.opt_tm.toInteger() > params.max_tm.toInteger() || params.opt_tm.toInteger() < params.min_tm.toInteger() ) {exit 1, " opt_tm: ${params.opt_tm} should > min_tm (${params.min_tm}) and < max_tm (${params.max_tm})"}
if (params.min_gc.toInteger() > params.max_gc.toInteger() ) {exit 1, " min_gc and max_gc: max_gc (${params.min_gc}) should be > min_gc(${params.max_gc})"}
if (params.opt_gc.toInteger() > params.max_gc.toInteger() || params.opt_gc.toInteger() < params.min_gc.toInteger() ) {exit 1, " opt_gc: ${params.opt_gc} should > min_gc (${params.min_gc}) and < max_gc (${params.max_gc})"}


// print run info

log.info """\
==============================================
CIRCprimerXL pipeline
==============================================
OncoRNALab - Marieke Vromman
https://github.com/OncoRNALab/CIRCprimerXL
https://hub.docker.com/repository/docker/oncornalab/CIRCprimerXL
==============================================
your input file : ${params.input_bed}
your output directory : ${params.output_dir}
"""

// define and run each process

process split_circRNAs {

	input:
	path(input_bed_handle)
	path(chrom_file_handle)

	output:
	path 'circ*'
	path 'start_time.txt'
	path 'all_circ.txt'

	"""
	01_validate_bed.py -i $input_bed_handle -c $chrom_file_handle
	02_split_circRNAs.py -i $input_bed_handle
	python3 -c 'from datetime import datetime; print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))' > start_time.txt
	"""

}


process get_seq {

	input:
	file(ind_circ_file_handle)
	path(fasta_index)
	path(exons)
	path(enst) 

	output:
	tuple val("${ind_circ_file_handle.baseName}"), path('input_primer3*')
	tuple val("${ind_circ_file_handle.baseName}"), path('input_NUPACK*')
	tuple val("${ind_circ_file_handle.baseName}"), path('input_SNP*')
	tuple val("${ind_circ_file_handle.baseName}"), path('input_filter*')
	tuple val("${ind_circ_file_handle.baseName}"), path('fasta_in.txt')
	tuple val("${ind_circ_file_handle.baseName}"), path('fasta_track.txt')
	tuple val("${ind_circ_file_handle.baseName}"), path('annotation*.txt')

	"""
	03_get_circ_seq_fastahack1.py -n $params.temp_l -i $ind_circ_file_handle -s $params.splice -e $exons -t $enst
	cat fasta_in.txt | /bin/fastahack-1.0.0/fastahack -c $fasta_index/$params.index_fasta_name > fasta_out.txt
	04_get_circ_seq_fastahack2.py -i $ind_circ_file_handle -n $params.primer3_diff -p $params.primer3_nr -a $params.min_tm -b $params.max_tm -c $params.opt_tm -d $params.diff_tm -e $params.min_gc -f $params.max_gc -g $params.opt_gc -j $params.amp_min -k $params.amp_max
	"""
}


process get_SNPs {

	container 'oncornalab/primerxl_circ:v0.11'

	input:
	tuple val(snp_id), path(in_SNP_handle)
	tuple val(snp_id), path(fasta_SNP_handle)
	tuple val(snp_id), path(fasta_track_handle)


	output:
	tuple val(snp_id), path('output_SNP*')
	tuple val(snp_id), path('output_SNP*')


	"""
	05_get_SNPs.py -i $in_SNP_handle -u $params.snp_url -f $fasta_SNP_handle -t $fasta_track_handle
	"""
	
}

process folding_template {

	input:
	tuple val(fold_id_t), path('in_folding_handle')

	output:
	tuple val (fold_id_t), path('output_NUPACK_*')
	tuple val (fold_id_t), path('output_NUPACK_*')

	"""
	06_get_sec_str_temp.py -i in_folding_handle
	"""
}



process get_primers {

	publishDir "$params.output_dir/primer3_details", mode: 'copy', pattern: 'output_primer3_*'
	errorStrategy "ignore"
	
	input:
	tuple val(u_filter_id), path('out_folding_template_upfront_filter_handle'), path('in_primer3_handle'), path('out_SNP_upfront_filter_handle')
	path(primer_settings_handle)
 
	output:
	tuple val(u_filter_id), path('amplicon_folding_input_circ*')
	tuple val(u_filter_id), path('all_primers_circ*')
	path('primer_spec_input_circ*')
	path('output_primer3_*')

	"""
	07_upfront_filter.py -i in_primer3_handle -a out_folding_template_upfront_filter_handle -f $params.upfront_filter -s out_SNP_upfront_filter_handle -l $params.temp_l
	/bin/primer3-2.5.0/src/primer3_core --output=output_primer3_${u_filter_id}.txt --p3_settings_file=$primer_settings_handle primer3_file*
	08_split_primers.py -i output_primer3_${u_filter_id}.txt
	"""
}

process folding_amplicon {

	input:
	tuple val(fold_id_a), path(amplicon_folding_in_handle)

	output:
	tuple val(fold_id_a), path('output_NUPACK_*')

	"""
	09_get_sec_str_amp.py -i $amplicon_folding_in_handle
	"""
}



process specificity_primers {

	container 'oncornalab/primerxl_circ:v0.11'

	// cpus 4
	// cpus 16

	input:
	path(bowtie_index)
	file('primer_spec_input_*')

	output:
	file('out_spec_primer.txt')

	script:
	"""
	cat primer_spec_input_* > primer_spec_input_all_circ.txt
	/bin/bowtie-1.3.0-linux-x86_64/bowtie --tryhard -X1000 -v3 --quiet -x $bowtie_index/$params.index_bowtie_name --threads ${task.cpus} --12 primer_spec_input_all_circ.txt > out_spec_primer.txt
	"""
}


process filter_primers {

	input:
	tuple val(filter_id), path('circ_file_handle'), path('out_SNP_handle'), path('out_folding_template_handle'), path('out_folding_amplicon_handle'), path('all_primers_per_circ_handle'), path('annotation_splice_handle')
	file(out_spec_primer_handle)
	
	output:
	path('selected_primers_*')
	path('all_primers')
	path('log_file_*')

	"""
	mkdir all_primers
	10_filter.py -A $circ_file_handle -P $all_primers_per_circ_handle -b $out_spec_primer_handle -l $params.temp_l -s $out_SNP_handle -t $out_folding_template_handle -a $out_folding_amplicon_handle -f $params.snp_filter -p $params.spec_filter
	11_gather_output.py -i all_primers/filtered_primers_* -a $annotation_splice_handle
	"""
}



process print_output {

	publishDir params.output_dir, mode: 'copy'

	input:
	file('results_per_circ*')
	file('log_file_per_circ*')
	file('start_time_file')
	path('all_primer_files')
	path('all_circ_file')


	output:
	path('all_primers')
	path('filtered_primers.txt')
	path('log_file.txt')
	path('summary_run.txt')

	"""
	mkdir all_primers
	cp all_primer_files*/* all_primers/
	echo "circ_ID	chr	start	end	primer_ID	FWD_primer	REV_primer	FWD_pos	FWD_length	REV_pos	REV_length	FWD_Tm	REV_Tm	FWD_GC	REV_GC	amplicon	PASS	start_annotation	end_annotation	splicing" > filtered_primers.txt
	cat results_per_circ* >> filtered_primers.txt
	echo "circ_ID	chr	start	end	design	primer_found	total_primer_pairs	passed	failed_spec	failed_SNP	failed_sec_str_temp	failed_sec_str_amp" > log_file.txt
	cat log_file_per_circ* >> log_file.txt
	12_summary_run.py -l log_file.txt -s start_time_file -o . -u $params.upfront_filter -a all_circ_file
	"""
}

workflow {
	split_circRNAs(
		input_bed,
		chrom_file)
	get_seq(
		split_circRNAs.out[0].flatten(),
		index_fasta,
		params.known_exons,
		params.list_ENST)
	get_SNPs(
		get_seq.out[2],
		get_seq.out[4],
		get_seq.out[5]
	)
	folding_template(
		get_seq.out[1]
	)
	get_primers(
		folding_template.out[0].join(get_seq.out[0]).join(get_SNPs.out[0]).groupTuple(),
		params.primer_settings
	)
	folding_amplicon (
		get_primers.out[0].transpose()
	)
	specificity_primers(
		index_bowtie,
		get_primers.out[2].collect()
	)
	filter_primers(
		get_seq.out[3].join(get_SNPs.out[1]).join(folding_template.out[1]).join(folding_amplicon.out[0]).join(get_primers.out[1]).join(get_seq.out[6]).groupTuple(),
		specificity_primers.out[0]
	)
	print_output(
		filter_primers.out[0].collect(),
		filter_primers.out[2].collect(),
		split_circRNAs.out[1],
		filter_primers.out[1].collect(),
		split_circRNAs.out[2]
	)

}