# CIRCprimerXL
Collaborators: Marieke Vromman, Pieter-Jan Volders

Questions concerning the GitHub structure/scripts can be addressed to any of the collaborators.

Primer design pipeline for circRNAs based on primerXL (Lefever, S., Pattyn, F., De Wilde, B. et al. High-throughput PCR assay design for targeted resequencing using primerXL. BMC Bioinformatics 18, 400 (2017). https://doi.org/10.1186/s12859-017-1809-3).

This pipeline runs entirly in the [oncornalab/primerxl_circ](https://hub.docker.com/repository/docker/oncornalab/primerxl_circ) docker image, which is available on DockerHub. It is not necessary to download this image locally, as Nextflow pulls the latest version automatically from DockerHub.

## Installation
### Required files and reference genomes
Two references are required to run the pipeline:
<ul>
  <li>a fasthack reference genome to extract the circRNA sequence surrounding the BSJ (https://github.com/ekg/fastahack)</li>
  <li>a Bowtie cDNA + ncRNA reference to test the specificity of the primers (https://github.com/BenLangmead/bowtie)</li>
    <li>a bed file containing all exons. This can either be downloaded using the UCSC Table Browser (https://genome.ucsc.edu/cgi-bin/hgTables?command=start, select group: all tables, table: knownGene, get output) or ensembl (https://uswest.ensembl.org/info/data/ftp/index.html, download the gtf file and transform into bed format using the gtf_to_bed.py python script (can be run in the Docker image)).</li>

</ul>

Both indexes below are examples and can be replaced by your choosen index (--index_fasta, --index_bowtie) and species (see below).

#### fastahack
Note: the fastahack index needs to be unzipped.

```bash
wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
fastahack -i Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

#### Bowtie
In general, a combination of the cDNA and ncRNA index are used to test the specificity of the primers.

```bash
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
gunzip Homo_sapiens.GRCh38.ncrna.fa.gz
cat Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.ncrna.fa > hg38_cdna.fa
rm Homo_sapiens.GRCh38.cdna.all.fa
rm Homo_sapiens.GRCh38.ncrna.fa
bowtie-build hg38_cdna.fa hg38_cdna
```

You do not need to install fastahack and Bowtie locally to create the required indexes. Instead you can use the Docker container. For this, [Docker](https://docs.docker.com/get-docker/) needs to be installed on your computer.
```bash
docker run -v "$PWD/assets":/assets oncornalab/primerxl_circ:v0.11 /bin/bowtie-1.3.0-linux-x86_64/bowtie-build /assets/index_bowtie/hg38_cdna.fa /assets/index_bowtie/hg38_cdna
docker run -v "$PWD/assets":/assets oncornalab/primerxl_circ:v0.11 /bin/fastahack-1.0.0/fastahack -i assets/index_fastahack/GRCh38_latest_genomic.fna
```


### Running on your computer
[Nextflow](https://www.nextflow.io/) and [Docker](https://docs.docker.com/get-docker/) should be installed locally. Make sure [Docker Desktop](https://www.docker.com/products/docker-desktop) is running when you want run the pipeline.

### Running on the HPC (UGent)
Nextflow version 20.10.0 is available on all clusters (swalot, skitty, victini, joltik, kirlia, doduo). The pipeline can be run through an interactive session. The pipeline can only run from the $VSC_SCRATCH_VO_USER directory.

```
qsub -I -l nodes=1:ppn=16 -l walltime=04:00:00
cd $VSC_SCRATCH_VO_USER/CIRCprimerXL/
module load Nextflow/20.10.0
nextflow run CIRCprimerXL.nf --help
```


## General usage

```
$ nextflow run CIRCprimerXL.nf --help

Usage:
	
	The typical command for running the pipeline is as follows:
	nextflow run CIRCprimerXL.nf -profile singularity
	
	Mandatory nextflow arguments:
	-profile            set to 'local' when running locally, set to 'singularity' when running on the HPC

	Mandatory pipeline arguments:
	--input_bed         path to input file with circRNAs in bed format (0-based annotation)
	--index_bowtie      path to bowtie genome index directory
	--index_bowtie_name the basename of the Bowtie index to be searched (the name of any of the index files up to but not including the final .1.ebwt / .rev.1.ebwt / ...)
	--index_fasta       path to fastahack index
	--index_fasta_name  the name of the fastahack genome file



	Optinal pipeline arguments:
	--splice			when set to 'yes' the input sequence will be spliced, when set to 'no' the input sequence will be unspliced
	--primer_settings   path to file with primer3plus settings (see primer3 manual)
	--chrom_file        file containing all chromosome sizes (to validate bed file)
	--known_exons       bed file containing exon annotation
	--list_ENST         file containing ENST numbers of canonical transcripts or transcripts of interest (this file can also be left empty)
	--primer3_diff      the minimum number of base pairs between the 3' ends of any two left primers (see also primer3 PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE)
	--primer3_nr        the number of primers designed by primer3; caution: setting this parameter to a large value will increase running time
	--min_tm	    minimum melt temperature of the primers (default: 58)
	--max_tm	    maximum melt temperature of the primers(default: 60)
	--opt_tm	    optimal melt temperature of the primers(default: 59)
	--diff_tm	    maximum difference in melt temperature between the primers(default: 2)
	--min_gc	    minimum GC contect of the primers (default: 30)
	--max_gc	    maximum GC contect of the primers(default: 80)
	--opt_gc	    optimal GC contect of the primers(default: 50)
	--amp_min	    minimum amplicon length (default: 60)
	--amp_max	    maximum amplicon length (default: 0)
	--temp_l            the number of nucleotides on each side of the circRNA BSJ that will be used for the template (example 150 => template of 300 nts in total)
	--spec_filter       when set to 'strict', only 2MM + 3MM are allowed; when set to 'loose', 2MM + 2MM and 2MM + 1MM are also allowed
	--snp_filter        when set to 'strict', no common SNPs are allowed in primer sequence; when set to 'loose', common SNPs are allowed in 5' first half of primer
	--snp_url           when using a differente species than human, the correct SNP database url should be provided; alternatively, this paramater can be set to 'off' if no SNP database is available
	--upfront_filter    when set to 'yes', SNPs and secundary structures are avoided before primer design; when set to 'str', secundary structures are avoided before primer design; when set to 'snp', snp are avoided before primer design; when set to 'no', no filtering before primer design is performed
	--output_dir        path to directory where the output files will be saved
	
```

This repository contains an example run with 3 circRNAs. For this, a small subset of the indexes are also present in the example folder. The example can be run by
```
nextflow CIRCprimerXL.nf -profile example
```
This is equivalent to
```
nextflow CIRCprimerXL.nf -profile local --output_dir example/output --input_bed example/input_circRNAs.bed --index_fasta example/index_fastahack --index_bowtie example/index_bowtie --index_bowtie_name = hg38_cdna_small
```

You can easily create your own profiles by modifying the nextflow.config file.

Nextflow keeps track of all the processes executed in your pipeline. If you want to restart the pipeline after a bug, add '-resume'. The execution of the processes that are not changed will be skipped and the cached result used instead.

Note: The pipeline results are cached by default in the directory $PWD/work. This folder can take of lot of disk space. If your are sure you won’t resume your pipeline execution, clean this folder periodically.

Note: If a circRNA is smaller than the requested template size, the template size is reduced to the circRNA size. Of note, if this 300-nucleotide template sequence includes an exon-intron boundary, the intronic region (which may not be part of the circRNA) is included. Some circRNAs effectively also include intronic sequences, and some BSJ concatenate an exonic and an intronic sequence. 

## Output
In the output folder, you will find
<ul>
  <li>filtered_primers.txt, a file containing one selected primer pair per circRNA (see below for column details)</li>
  <li>log_file.txt, a file </li>
  <li>summary_run.txt, </li>
  <li>all_primers directory</li>
  <li>primer3_details directory</li>
</ul>

filtered_primers.txt output file column names:

| column name      | description                                                                                                            |
|:-----------------|:-----------------------------------------------------------------------------------------------------------------------|
| circ_ID          | circ id assigned to each circRNA (unique within one run)                                                               |
| chr              | circRNA chromome                                                                                                       |
| start            | circRNA start position                                                                                                 |
| end              | circRNA end position                                                                                                   |
| primer_ID        | primer ID generated by primer3                                                                                         |
| FWD_primer       | forward primer                                                                                                         |
| REV_primer       | reverse primer                                                                                                         |
| FWD_pos          | relative position of forward primer                                                                                    |
| FWD_length       | length of forward primer                                                                                               |
| REV_pos          | relative position of reverse primer                                                                                    |
| REV_length       | length of reverse primer                                                                                               |
| FWD_Tm           | melt temperature of forward primer                                                                                     |
| REV_Tm           | melt temperature of reverse primer                                                                                     |
| FWD_GC           | GC content of forward primer                                                                                           |
| REV_GC           | GC content of reverse primer                                                                                           |
| amplicon         | amplicon sequence amplified by the primer pair                                                                         |
| PASS             | result of filtering (PASS if the primer pair passed all filters, FAIL if   the primer pair failed one or more filters) |
| start_annotation | exons used for the template sequence on the right side of the BSJ                                                      |
| end_annotation   | exons used for the template sequence on the left side of the BSJ                                                       |
| splicing         | spliced or unspliced template sequence was used                                                                        |


## Other species
As default, CIRCprimerXL designs primers voor humans. To design primers for other species, the following files have to be provided and parsed through the corresponding parameters:
<ul>
  <li>a file containing the chromosome sizes (parameter chrom_file) (for example from: https://www.ncbi.nlm.nih.gov/grc/human/data)</li>
  <li>a fastahack index (parameter index_fasta) and Bowtie index (parameter index_bowtie) (see also above)</li>
  <li>a SNP database link (parameter snp_url) (for example: http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153Common.bb
)</li>
  <li>a file containing all ENST numbers of canonical transcripts or transcrits. This file can be generated by downloading the MANE file (https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/) and transforming it into a simple list using generate_ENST_list.py (can be run in the Docker image).
</ul>


## Nextflow tower

[Nextflow tower](https://tower.nf/) can be used to monitor the pipeline while it's running.
```
nextflow run CIRCprimerXL.nf -with-tower
```

When Nextflow tower is used in combination with the HPC, the nextflow version and tower access token should be indicated.
```
export NXF_VER="20.10.0"
export TOWER_ACCESS_TOKEN=your_token_here
```

