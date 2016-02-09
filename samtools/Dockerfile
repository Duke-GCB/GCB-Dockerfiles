FROM phusion/baseimage
MAINTAINER Alejandro Barrera <alejandro.barrera@duke.edu>

# Install dependencies
RUN apt-get update && apt-get install -y \
  curl \
  build-essential \
  libncurses-dev \ 
  zlib1g-dev

# Installs samtools from compiled distribution into /opt/bowtie
ENV SAMTOOLS_RELEASE=1.3
ENV SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/
ENV DEST_DIR=/opt/

# Download Samtools; untar & decompress; remove tr.bz2 file; compile & install samtools; remove unnecessary files
RUN curl -SLo ${DEST_DIR}/samtools-${SAMTOOLS_RELEASE}.tar.bz2 ${SAMTOOLS_URL}/${SAMTOOLS_RELEASE}/samtools-${SAMTOOLS_RELEASE}.tar.bz2 && \
  tar -xf ${DEST_DIR}/samtools-${SAMTOOLS_RELEASE}.tar.bz2 -C ${DEST_DIR} && \
  rm ${DEST_DIR}/samtools-${SAMTOOLS_RELEASE}.tar.bz2 && \
  cd ${DEST_DIR}/samtools-${SAMTOOLS_RELEASE} && \
  ./configure && \
  make && \
  make install && \
  rm -rf ${DEST_DIR}/samtools-${SAMTOOLS_RELEASE}

CMD ["samtools"]
