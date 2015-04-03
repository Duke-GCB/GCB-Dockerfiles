FROM phusion/baseimage
MAINTAINER Dan Leehr <dan.leehr@duke.edu>

# Install fastqc from apt

RUN apt-get update && apt-get install -y \
  openjdk-7-jre-headless \
  fastqc

# Copy wrapper script

COPY fastqc.sh /usr/bin/fastqc.sh

CMD ["/usr/bin/fastqc.sh"]