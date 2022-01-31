# based on existing java image
FROM oncornalab/primerxl_circ:v0.26


# copied from https://hub.docker.com/r/picoded/ubuntu-openjdk-8-jdk/dockerfile/
RUN apt-get update -y --fix-missing
# This is in accordance to : https://www.digitalocean.com/community/tutorials/how-to-install-java-with-apt-get-on-ubuntu-16-04
RUN apt-get update && \
	apt-get install -y openjdk-8-jdk && \
	apt-get install -y ant && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/* && \
	rm -rf /var/cache/oracle-jdk8-installer;
	
# Fix certificate issues, found as of 
# https://bugs.launchpad.net/ubuntu/+source/ca-certificates-java/+bug/983302
RUN apt-get update && \
	apt-get install -y ca-certificates-java && \
	apt-get clean && \
	update-ca-certificates -f && \
	rm -rf /var/lib/apt/lists/* && \
	rm -rf /var/cache/oracle-jdk8-installer;

# Setup JAVA_HOME, this is useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME

# get nextflow
WORKDIR /usr/bin
RUN wget -qO- https://get.nextflow.io | bash
RUN chmod +x nextflow
RUN nextflow run hello -plugins nf-amazon@1.3.4

# get pipeline
ADD . /CIRCprimerXL
WORKDIR /CIRCprimerXL
ENTRYPOINT ["nextflow", "/CIRCprimerXL/CIRCprimerXL.nf", "-config", "/CIRCprimerXL/nf-runner.config"]
