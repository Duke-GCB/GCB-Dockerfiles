#!/bin/bash

set -e

# This is a script for a docker image to run trimmomatic, a java program
# See http://www.usadellab.org/cms/?page=trimmomatic
# It translates environment variables and wraps fastqc commands
# Paired mode
#
# CONT prefix means this variable is expected to be set in the container at runtime

# INPUTS
#
# CONT_INPUT_FASTQ_READ_1
# CONT_INPUT_FASTQ_READ_2
# CONT_INPUT_ADAPTERS_FILE
#
# PARAMETERS
#
# CONT_PARAM_QUALITY_SCORE  default -phred33
# CONT_PARAM_ADAPTERS_PARAMS
# CONT_PARAM_STEPS default  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30
# OUTPUTS
#
# CONT_OUTPUT_PAIRED_1
# CONT_OUTPUT_UNPAIRED_1
# CONT_OUTPUT_PAIRED_2
# CONT_OUTPUT_UNPAIRED_2
# CONT_OUTPUT_ERROR_LOG

# PE is paired ends mode, used for RNA_SEQ
# command from IMD
# java -jar /data/reddylab/software/Trimmomatic-0.32/trimmomatic-0.32.jar PE \
#   -threads 1 \                # threads
#   -phred33 \                  # quality score - Phred-33
#   RNA_seq_1.fastq \           # input1
#   RNA_seq_2.fastq \           # input2
#   RNA_seq_1.trimmed.P.fastq \ # paired output 1
#   RNA_seq_1.trimmed.U.fastq \ # unpaired output 1
#   RNA_seq_2.trimmed.P.fastq \ # paired output 2
#   RNA_seq_2.trimmed.U.fastq \ # unpaired output 2
#   adapters.fa:2:30:15 \       # Are we missing ILLUMINACLIP here?
#   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30 \ # more steps
#   2 > RNA.trim_err.txt # Redirect stderr

# old code below here

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
