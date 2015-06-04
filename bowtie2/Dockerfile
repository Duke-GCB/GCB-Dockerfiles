FROM phusion/baseimage
MAINTAINER Dan Leehr <dan.leehr@duke.edu>

# Install dependencies

RUN apt-get update && apt-get install -y \
  bowtie2

# Copy wrapper script

COPY bowtie2.sh /usr/bin/bowtie2.sh

CMD ["/usr/bin/bowtie2.sh"]
