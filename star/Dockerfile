FROM phusion/baseimage
MAINTAINER Dan Leehr <dan.leehr@duke.edu>

# Install git

RUN apt-get update && apt-get install -y \
  git \
  build-essential \
  curl

# Copy wrapper script

COPY star.sh /usr/bin/star.sh

CMD ["/usr/bin/star.sh"]
