FROM ubuntu:16.04
LABEL maintainer="dan.leehr@duke.edu"

# SeekDeep installation based on http://baileylab.umassmed.edu/SeekDeep/installingSeekDeep_Ubuntu.html#brand-new-ubuntu-machine

# RUN echo America/New_York | tee  /etc/timezone  && dpkg-reconfigure --frontend noninteractive tzdata
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
  build-essential \
  software-properties-common \
  python-software-properties \
  curl \
  git \
  make
RUN add-apt-repository ppa:ubuntu-toolchain-r/test -y
RUN apt-get update && apt-get install -y g++-7

# Ubuntu base image not configured with a locale. Without this, the seekdeep build fails
# with cryptic python UnicodeDecodeErrors (trying to decode to a strange default locale)
# Via https://stackoverflow.com/a/38553499

RUN apt-get install -y locales
RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    dpkg-reconfigure --frontend=noninteractive locales && \
    update-locale LANG=en_US.UTF-8
ENV LANG=en_US.UTF-8

# SeekDeep requires CMake 3.7
RUN apt-get install -y cmake

ENV SEEKDEEP_RELEASE=2.5.2
ENV SEEKDEEP_URL=https://github.com/bailey-lab/SeekDeep/archive/
ENV DEST_DIR=/opt/SeekDeep

RUN mkdir ${DEST_DIR}
WORKDIR ${DEST_DIR}

# Download and extract
RUN curl -SLo ${DEST_DIR}/SeekDeep-${SEEKDEEP_RELEASE}.tar.gz ${SEEKDEEP_URL}/v${SEEKDEEP_RELEASE}.tar.gz && \
  tar -xzf ${DEST_DIR}/SeekDeep-${SEEKDEEP_RELEASE}.tar.gz -C ${DEST_DIR} && \
  rm ${DEST_DIR}/SeekDeep-${SEEKDEEP_RELEASE}.tar.gz

# Install
WORKDIR ${DEST_DIR}/SeekDeep-${SEEKDEEP_RELEASE}
RUN ./setup.py --symlinkBin
RUN ./setup.py --addBashCompletion
RUN ./install.sh 2

# Put seekDeep in path
ENV PATH ${DEST_DIR}/SeekDeep-${SEEKDEEP_RELEASE}/bin:$PATH

CMD ["SeekDeep"]
