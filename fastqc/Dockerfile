FROM ubuntu:12.04
MAINTAINER Dan Leehr <dan.leehr@duke.edu>

# fastqc requires java
RUN apt-get update && apt-get install -y \
  curl \
  unzip \
  perl \
  openjdk-7-jre-headless

# Installs fastqc from compiled java distribution into /opt/FastQC
ENV FASTQC_URL http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
ENV FASTQC_RELEASE fastqc_v0.11.4.zip
ENV DEST_DIR /opt/

# Make destination directory
RUN mkdir -p $DEST_DIR

# Download & extract archive - Repo includes binaries for linux
WORKDIR /tmp

# Do this in one command to avoid caching the zip file and its removal in separate layers
RUN curl -SLO ${FASTQC_URL}/${FASTQC_RELEASE} && unzip ${FASTQC_RELEASE} -d ${DEST_DIR} && rm ${FASTQC_RELEASE}

# Make the wrapper script executable
RUN chmod +x ${DEST_DIR}/FastQC/fastqc

# Include it in PATH
ENV PATH ${DEST_DIR}/FastQC:$PATH
