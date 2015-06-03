#!/bin/bash

set -e

# This is a script for a docker image to run star to generate a genome index

# INPUTS
#
# CONT_INPUT_SJDB_GTF_FILE
# CONT_INPUT_GENOME_FASTA_FILES

# OUTPUTS
#
# CONT_OUTPUT_DIR - output directory for temp files and logs
# CONT_OUTPUT_GENOME_DIR - directory for generating STAR-formatted genome

# PARAMS
#
# CONT_PARAM_THREADS, default 32
# CONT_PARAM_SJDB_OVERHANG, default 25

# Check that variables are set
[ -z "$CONT_INPUT_SJDB_GTF_FILE" ] && echo "Error: The CONT_INPUT_SJDB_GTF_FILE variable must be set" && exit 1
[ -z "$CONT_INPUT_GENOME_FASTA_FILES" ] && echo "Error: The CONT_INPUT_GENOME_FASTA_FILES variable must be set" && exit 1
[ -z "$CONT_OUTPUT_DIR" ] && echo "Error: The CONT_OUTPUT_DIR variable must be set" && exit 1
[ -z "$CONT_OUTPUT_GENOME_DIR" ] && echo "Error: The CONT_OUTPUT_GENOME_DIR variable must be set" && exit 1

# Check that output directories are writable
[ ! -w "$CONT_OUTPUT_DIR" ] && echo "Error: output dir $CONT_OUTPUT_DIR is not writable" && exit 1
[ ! -w "$CONT_OUTPUT_GENOME_DIR" ] && echo "Error: output dir $CONT_OUTPUT_GENOME_DIR is not writable" && exit 1

# Populate defaults
THREADS="32"
SJDB_OVERHANG="25"

if [ ! -z "$CONT_PARAM_THREADS" ]; then
  THREADS=$CONT_PARAM_THREADS
fi

if [ ! -z "$CONT_PARAM_SJDB_OVERHANG" ]; then
  SJDB_OVERHANG=$CONT_PARAM_SJDB_OVERHANG
fi

# Build command
STAR_BIN=$(which STAR)

STAR_CMD="$STAR_BIN \
  --outFileNamePrefix $CONT_OUTPUT_DIR \
  --runMode genomeGenerate \
  --sjdbGTFfile $CONT_INPUT_SJDB_GTF_FILE \
  --sjdbOverhang $SJDB_OVERHANG \
  --genomeDir $CONT_OUTPUT_GENOME_DIR \
  --genomeFastaFiles $CONT_INPUT_GENOME_FASTA_FILES \
  --runThreadN $THREADS"

echo
echo "Starting $0..."
echo "$STAR_VERSION"
echo "Running star:"
echo "$STAR_CMD"
sh -c "$STAR_CMD"
