#!/bin/bash

# CI Build Script for GCB-Dockerfiles
#
# Builds docker images when a commit to the repo changes a Dockerfile
# Tags versioned images as `latest` when a symlink `latest` points to the version's directory
#
# Detecting changes to build is done using git diff and comparing filenames.
# The build's checked-out branch causes the comparison to happen in one of two ways:
#
# 1. When building (and deploying) the deploy branch (typically 'master'), the script
#    expects the HEAD commit to be a merge commit (from merging a pull request).
#    If it is not a merge, the build fails.
#    Otherwise, the parents of the merge are compared to get the changed files.
#    Using this approach, we test the changes that are merged into the deploy branch.
#
# 2. When building other branches, the script compares the branch with the deploy
#    branch (e.g. 'master'). Any files changed between the two branches are considered
#    for building. This way, we can trigger rebuilding a docker image without adding
#    irrelevant comments or blank lines to Dockerfiles.

set -e
source functions.sh
# Check that the Docker Hub organization to use is in the DOCKERHUB_ORG variable
check_org
# Check that the branch to use for deploying (typically 'master') is in the DEPLOY_BRANCH variable
check_deploy_branch
# Get the range of commits to compare for detecting changed files
compare_range=$(get_compare_range)
# Get a list of changed paths in the repo to look for <tool>/<version>/Dockerfile
paths=$(changed_paths_in_range "$compare_range")
# Print out what changes are being considered
print_changed "$compare_range" "$paths"
# Loop through the changed files and build Docker images for any that match
# <tool>/<version>/Dockerfile. If none found, prints a message indicating so.
build_images "$DOCKERHUB_ORG" "$paths"
