#!/bin/bash

set -e

# This is a script for a docker image to run fastqc, a java program
# See http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# It translates environment variables and wraps fastqc commands
#
# CONT prefix means this variable is expected to be set in the container at runtime

# INPUTS
#
# CONT_INPUT_FASTQ_FILE - a FASTQ file to process

# OUTPUTS
#
# CONT_OUTPUT_FASTQC_REPORT - The FastQC report files, in zip format

# Check that variables are set
[ -z "$CONT_INPUT_FASTQ_FILE" ] && echo "Error: The CONT_INPUT_FASTQ_FILE variable must be set" && exit 1

# Check that output file does not exist
# These break if there are spaces
[ -e "$CONT_OUTPUT_FASTQC_REPORT" ] && echo "Error: output file at $CONT_OUTPUT_FASTQC_REPORT already exists" && exit 1

# Check that output file is writable
[ ! -w $(dirname "$CONT_OUTPUT_FASTQC_REPORT") ] && echo "Error: output file $CONT_OUTPUT_FASTQC_REPORT is not writable" && exit 1

FASTQC_BIN=$(which fastqc)
FASTQC_VERSION=$($FASTQC_BIN --version)
WORKDIR=$(mktemp -d)
FASTQC_CMD="$FASTQC_BIN -o $WORKDIR $CONT_INPUT_FASTQ_FILE --noextract"

echo
echo "Starting $0..."
echo "$FASTQC_VERSION"
echo "$FASTQC_CMD"
sh -c "$FASTQC_CMD"

# Output files are created in workdir, and are named basename of input file + _fastqc.zip
INPUT_FILE_NAME=$(basename $CONT_INPUT_FASTQ_FILE)
ZIPFILE=$WORKDIR/${INPUT_FILE_NAME%.*}_fastqc.zip
mv $ZIPFILE $CONT_OUTPUT_FASTQC_REPORT
