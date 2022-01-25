# based on existing Docker image
FROM ubuntu:20.04

# tell Docker it's not interactive
ARG DEBIAN_FRONTEND=noninteractive

# install 'make' function (in build-essential & others from primer3 manual)
RUN apt-get update -y --fix-missing
RUN apt-get install -y --fix-missing build-essential g++ cmake git-all
RUN apt-get install -y wget
RUN apt-get install -y unzip
RUN apt-get install -y python3
RUN apt-get install -y python3-pip
RUN apt-get install -y bedtools

# install python packages
RUN python3 -m pip install -U pip Bio
RUN python3 -m pip install -U matplotlib jupyterlab
RUN python3 -m pip install -U pandas
RUN python3 -m pip install -U datetime
RUN python3 -m pip install -U pybedtools


# install NUPACK
ADD assets/nupack-4.0.0.23 /bin/nupack-4.0.0.23

WORKDIR /bin

RUN python3 -m pip install -U nupack -f ./nupack-4.0.0.23/package
RUN python3 -m pip install -U pytest
RUN python3 -m pytest -v


# install primer3 version 2.5.0
RUN wget https://github.com/primer3-org/primer3/archive/v2.5.0.tar.gz
RUN tar -xvf v2.5.0.tar.gz
WORKDIR primer3-2.5.0/src
RUN make all

WORKDIR /

# install dependecies of bigBedtoBed
RUN apt-get install -y libkrb5-3
RUN apt-get install -y curl

# install bigBedtoBed
WORKDIR /bin
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
RUN chmod a+x bigBedToBed

WORKDIR /

# install bowtie
WORKDIR /bin
RUN wget https://github.com/BenLangmead/bowtie/releases/download/v1.3.0/bowtie-1.3.0-linux-x86_64.zip
RUN unzip bowtie-1.3.0-linux-x86_64.zip
WORKDIR /

# install fastahack
WORKDIR /bin
RUN wget https://github.com/ekg/fastahack/archive/refs/tags/v1.0.0.zip
RUN unzip v1.0.0.zip
WORKDIR fastahack-1.0.0
RUN make all
WORKDIR /

# install python scripts
ADD scripts /usr/bin
RUN chmod +x /usr/bin/get_sec_str_temp.py /usr/bin/get_sec_str_amp.py /usr/bin/filter.py /usr/bin/split_primers.py /usr/bin/split_circRNAs.py /usr/bin/gather_output.py /usr/bin/get_SNPs.py /usr/bin/upfront_filter.py /usr/bin/validate_bed.py /usr/bin/summary_run.py /usr/bin/get_circ_seq_fastahack2.py usr/bin/get_circ_seq_fastahack1.py /usr/bin/gtf_to_bed.py /usr/bin/generate_ENST_list.py
