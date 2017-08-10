FROM ubuntu:14.04
MAINTAINER Dan Leehr <dan.leehr@duke.edu>

# Install dependencies
ENV BOWTIE2_RELEASE=2.1.0-2

RUN apt-get update && apt-get install -y \
  bowtie2=${BOWTIE2_RELEASE}

CMD ["bowtie2"]
