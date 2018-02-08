#!/bin/bash

# CI Deploy Script for GCB-Dockerfiles
#
# Expects to be authenticated to Docker Hub and only run from master branch
# Pushes Docker images and tags to Docker Hub when build-ci.sh has built the image

set -e
source functions.sh
# See build-ci.sh for explanation of these conventions/rules
check_org
sha=$(last_nonmerge_commit_sha)
paths=$(changed_paths_in_commit $sha)
push_images "$DOCKERHUB_ORG" "$paths"

