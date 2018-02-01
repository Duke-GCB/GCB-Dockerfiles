#!/bin/bash

# CI Build Script for GCB-Dockerfiles
#
# Builds docker images when a commit to the repo changes a Dockerfile
# Tags versioned images as `latest` when a symlink `latest` points to the version's directory

set -e
source functions.sh
# Docker Hub organization to prefix on built docker images should be in DOCKERHUB_ORG
check_org
# We only want to build images when a Dockerfile changes, so we start with a list of
# changed paths and look for <tool>/<version>/Dockerfile
# To get the list of changed paths, we use `git diff-tree`, which returns nothing for merge commits
# So we find the most recent commit that's not a merge.
sha=$(last_nonmerge_commit_sha)
# Get a list of changed paths in the repo to look for <tool>/<version>/Dockerfile
paths=$(changed_paths_in_commit $sha)
# Print out what changes are being considered
print_changed "$sha" "$paths"
# Loop through the changed files and build Docker images for any that match
# <tool>/<version>/Dockerfile. If none found, prints a message indicating so.
build_images "$DOCKERHUB_ORG" "$paths"
