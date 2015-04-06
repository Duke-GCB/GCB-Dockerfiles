FROM phusion/baseimage
MAINTAINER Dan Leehr <dan.leehr@duke.edu>

RUN apt-get update && apt-get install -y \
  openjdk-7-jre-headless \
  trimmomatic

# Copy wrapper script

COPY trimmomatic.sh /usr/bin/trimmomatic.sh

CMD ["/usr/bin/trimmomatic.sh"]