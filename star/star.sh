#!/bin/bash

set -e

# This is a script for a docker image to run star to generate a genome index

# INPUTS
#
# CONT_INPUT_SJDB_GTF_FILE
# CONT_INPUT_GENOME_FASTA_FILES

# OUTPUTS
#
# CONT_OUTPUT_GENOME_DIR

# PARAMS
#
# CONT_PARAM_THREADS, default 4
# CONT_PARAM_SJDB_OVERHANG, default 25

# Check that variables are set
[ -z "$CONT_INPUT_SJDB_GTF_FILE" ] && echo "Error: The CONT_INPUT_SJDB_GTF_FILE variable must be set" && exit 1
[ -z "$CONT_INPUT_GENOME_FASTA_FILES" ] && echo "Error: The CONT_INPUT_SJDB_GTF_FILE variable must be set" && exit 1
[ -z "$CONT_OUTPUT_GENOME_DIR" ] && echo "Error: The CONT_OUTPUT_GENOME_DIR variable must be set" && exit 1

# Check that output directory is writable
[ ! -w $ "$CONT_OUTPUT_GENOME_DIR") ] && echo "Error: output dir $CONT_OUTPUT_GENOME_DIR is not writable" && exit 1

# Populate defaults
THREADS="4"
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
  --runMode genomeGenerate \
  --sjdbGTFfile $CONT_INPUT_SJDB_GTF_FILE \
  --sjdbOverhang $SJDB_OVERHANG \
  --genomeDir $CONT_OUTPUT_GENOME_DIR \
  --genomeFastaFiles $CONT_INPUT_GENOME_FASTA_FILES \
  --runThreadN $THREADS

echo
echo "Starting $0..."
echo "$STAR_VERSION"
echo "Running star:"
echo "$STAR_CMD"
sh -c "$STAR_CMD"
