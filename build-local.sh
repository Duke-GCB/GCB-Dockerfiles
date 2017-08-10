#!/bin/bash

# Script for local development building of Docker images

ORG="dukegcb"

if [ "$1" == "" ]; then
  echo "Builds a docker image locally"
  echo "Usage: $0 <tool_dir>/<version_dir>"
  exit 1
fi

TOOL_VERSION=$1
TOOL=$(echo $TOOL_VERSION | cut -d "/" -f 1)
VERSION=$(echo $TOOL_VERSION | cut -d "/" -f 2)

docker build -f ${TOOL_VERSION}/Dockerfile -t ${ORG}/${TOOL}:${VERSION} ${TOOL_VERSION}
