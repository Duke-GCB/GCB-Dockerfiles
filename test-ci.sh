#!/bin/bash

set -e
source functions.sh
sha=$(last_nonmerge_commit_sha)
paths=$(changed_paths_in_commit $sha)
python3 tests/imagecheck.py "$DOCKERHUB_ORG" $paths
