FROM phusion/baseimage
MAINTAINER Dan Leehr <dan.leehr@duke.edu>

# Install trimmomatic from ...
# TODO: Where to get trimmomatic

RUN apt-get update && apt-get install -y \
  openjdk-7-jre-headless

# Copy wrapper script

COPY trimmomatic.sh /usr/bin/trimmomatic.sh

CMD ["/usr/bin/trimmomatic.sh"]