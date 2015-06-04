#!/bin/bash

set -e

# This is a script for a docker image to run bowtie2
# It translates environment variables and wraps bowtie2 commands
#
# CONT prefix means this variable is expected to be set in the container at runtime

# INPUTS
#
# CONT_INPUT_GENOME_INDEX_PREFIX - the prefix to the bt2 index file (minus trailing .X.bt2).
# CONT_INPUT_READ_FILE

# PARAMETERS
#
# CONT_PARAM_THREADS, default 8

# OUTPUTS
#
# CONT_OUTPUT_ALIGNED_FILE
# CONT_OUTPUT_ALIGNMENT_METRICS_FILE (optional)

# Check that variables are set
[ -z "$CONT_INPUT_GENOME_INDEX_PREFIX" ] && echo "Error: The CONT_INPUT_GENOME_INDEX_PREFIX variable must be set" && exit 1
[ -z "$CONT_INPUT_READ_FILE" ] && echo "Error: The CONT_INPUT_READ_FILE variable must be set" && exit 1
[ -z "$CONT_OUTPUT_ALIGNED_FILE" ] && echo "Error: CONT_OUTPUT_ALIGNED_FILE variable must be set" && exit 1

# Check that output file does not exist
[ -e "$CONT_OUTPUT_ALIGNED_FILE" ] && echo "Error: output file at $CONT_OUTPUT_ALIGNED_FILE already exists" && exit 1

# Check that output file is writable
OUTPUT_DIR=$(dirname "$CONT_OUTPUT_ALIGNED_FILE")
[ ! -w "$OUTPUT_DIR" ] && echo "Error: output file $CONT_OUTPUT_ALIGNED_FILE is not writable" && exit 1

# Populate defaults
THREADS="8"
METRICS="/dev/null"

if [ ! -z "$CONT_PARAM_THREADS" ]; then
  THREADS=$CONT_PARAM_THREADS
fi

if [ ! -z "$CONT_OUTPUT_ALIGNMENT_METRICS_FILE" ]; then
  METRICS=$CONT_OUTPUT_ALIGNMENT_METRICS_FILE
fi
# Build command
BOWTIE2_BIN=$(which bowtie2)
# Get version
BOWTIE2_VERSION=$($BOWTIE2_BIN --version)
BOWTIE2_CMD="$BOWTIE2_BIN \
  --threads $THREADS \
  -x $CONT_INPUT_GENOME_INDEX_PREFIX \
  -U $CONT_INPUT_READ_FILE \
  -S $CONT_OUTPUT_ALIGNED_FILE \
  2> $METRICS"

echo
echo "Starting $0..."
echo "$BOWTIE2_VERSION"
echo "Running bowtie2:"
echo "$BOWTIE2_CMD"
sh -c "$BOWTIE2_CMD"
