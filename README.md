GGB-Dockerfiles
==========

Docker images for bioinformatic workflows that are built on [dockerhub](https://hub.docker.com/) under the [dukegcb organization](https://hub.docker.com/u/dukegcb/).

## Source Repository Organization

Each top level directory in the [GCB-Dockerfiles](https://github.com/Duke-GCB/GCB-Dockerfiles/) repository
corresponds to a `dukegcb/<directory-name>` docker container. 

Within the top level directory each Dockerfile will be named like so: `Dockerfile-<software-version>`. 

For example a Dockerfile that builds version 0.11.4 of FastQC would be stored in a file at `fastqc/Dockerfile-0.11.4`.

## Docker Hub Setup
Each top level directory will have a Automated Build setup with dockerhub.

Under `Build Setttings` There will be a build rule latest and each version of the software pointing to the specific `<directory-name>Dockerfile-<software-version>` Dockefile on the master branch with the `<software-version>` as the tag name.

For example if we want to build version 0.11.4 and 0.11.5 of FastQC there would be two dockerfiles: `fastqc/Dockerfile-0.11.4` and `fastqc/Dockerfile-0.11.5`. 

We would create three build rules on dockerhub for the dukegcb/fastqc container:
```
Type    Name   Dockerfile Location         Docker Tag Name
Branch  master /fastqc/Dockerfile-0.11.5   latest
Branch  master /fastqc/Dockerfile-0.11.5   0.11.5
Branch  master /fastqc/Dockerfile-0.11.4   0.11.4
```
