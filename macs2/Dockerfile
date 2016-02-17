FROM phusion/baseimage
MAINTAINER Alejandro Barrera <alejandro.barrera@duke.edu>

# Install python (+devs) and essentials for building
RUN apt-get update && apt-get install -y \
        build-essential \
        python \
        python-numpy \
        python-dev \
        cython \
        python-setuptools && \
    apt-get clean autoclean && \
    apt-get autoremove -y

ENV MACS2_VERSION=2.1.0.20151222
ENV MACS2_URL=https://pypi.python.org/packages/source/M/MACS2/MACS2-${MACS2_VERSION}.tar.gz
ENV DEST_DIR=/opt/

# Download Macs, decompress, install python setup.py and clean unnecessary files
RUN curl -SLo ${DEST_DIR}/macs2-${MACS2_VERSION}.tar.gz ${MACS2_URL} && \
    tar -xf ${DEST_DIR}/macs2-${MACS2_VERSION}.tar.gz -C ${DEST_DIR} && \ 
    rm ${DEST_DIR}/macs2-${MACS2_VERSION}.tar.gz && \
    cd ${DEST_DIR}/MACS2-${MACS2_VERSION} && \
    python setup.py install && \
    rm -rf ${DEST_DIR}/MACS2-${MACS2_VERSION}

ENV PYTHON_EGG_CACHE=/tmp/
# Default command to execute at startup of the container
CMD ["macs2"]
