FROM r-base:latest

RUN apt-get update && apt-get install -y \
  libboost-all-dev \
  curl \
  build-essential \
  libncurses-dev \
  zlib1g-dev \
  gawk

# Installs samtools from sources into /opt/
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

COPY . /usr/local/src/myscripts
WORKDIR /usr/local/src/myscripts
RUN R CMD INSTALL bitops_1.0-6.tar.gz
RUN R CMD INSTALL caTools_1.17.1.tar.gz
RUN R CMD INSTALL spp_1.10.1.tar.gz
CMD ["Rscript", "run_spp.R"]