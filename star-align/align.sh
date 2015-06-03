#!/bin/bash

set -e

# This is a script for a docker image to run star to align

# INPUTS
#
# CONT_INPUT_GENOME_DIR
# CONT_INPUT_READ_FILE_1
# CONT_INPUT_READ_FILE_2 [optional]

# OUTPUTS
#
# CONT_OUTPUT_WORKDIR
# CONT_OUTPUT_ALIGNED_FILE - The aligned file

# PARAMS
#
# CONT_PARAM_THREADS, default 4
# CONT_PARAM_GENOME_LOAD, default LoadAndKeep
# CONT_PARAM_SAM_ATTRIBUTES, default NH HI AS NM MD
# CONT_PARAM_FILTER_INTRON_MOTIFS, default RemoveNonCanonical

# Need to add rest of options

# Check that variables are set
[ -z "$CONT_INPUT_GENOME_DIR" ] && echo "Error: The CONT_INPUT_GENOME_DIR variable must be set" && exit 1
[ -z "$CONT_INPUT_READ_FILE_1" ] && echo "Error: The CONT_INPUT_READ_FILE_1 variable must be set" && exit 1
[ -z "$CONT_OUTPUT_WORKDIR" ] && echo "Error: The CONT_OUTPUT_WORKDIR variable must be set" && exit 1
[ -z "$CONT_OUTPUT_ALIGNED_FILE" ] && echo "Error: CONT_OUTPUT_ALIGNED_FILE variable must be set" && exit 1

# Check that output file does not exist
[ -e "$CONT_OUTPUT_ALIGNED_FILE" ] && echo "Error: output file at $CONT_OUTPUT_ALIGNED_FILE already exists" && exit 1

# Check that output directory is writable
[ ! -w "$CONT_OUTPUT_WORKDIR" ] && echo "Error: working dir $CONT_OUTPUT_WORKDIR is not writable" && exit 1
# Check that output file is writable
[ ! -w $(dirname "$CONT_OUTPUT_ALIGNED_FILE") ] && echo "Error: output file $CONT_OUTPUT_ALIGNED_FILE is not writable" && exit 1

# Populate defaults
THREADS="4"
SAM_ATTRIBUTES="NH HI AS NM MD"
FILTER_INTRON_MOTIFS="RemoveNoncanonical"

if [ ! -z "$CONT_PARAM_THREADS" ]; then
  THREADS=$CONT_PARAM_THREADS
fi

if [ ! -z "$CONT_PARAM_SAM_ATTRIBUTES" ]; then
  SAM_ATTRIBUTES=$CONT_PARAM_SAM_ATTRIBUTES
fi

if [ ! -z "$CONT_PARAM_FILTER_INTRON_MOTIFS" ]; then
  FILTER_INTRON_MOTIFS=$CONT_PARAM_FILTER_INTRON_MOTIFS
fi

# Build command
STAR_BIN=$(which STAR)

STAR_CMD="$STAR_BIN \
  --outFileNamePrefix $CONT_OUTPUT_WORKDIR/ \
  --runThreadN $THREADS \
  --genomeDir $CONT_INPUT_GENOME_DIR \
  --readFilesIn $CONT_INPUT_READ_FILE_1 $CONT_INPUT_READ_FILE_2 \
  --outSAMattributes $SAM_ATTRIBUTES \
  --outFilterIntronMotifs $FILTER_INTRON_MOTIFS"

echo
echo "Starting $0..."
echo "$STAR_VERSION"
echo "Running star:"
echo "$STAR_CMD"
sh -c "$STAR_CMD"

# Generates alignment as $CONT_OUTPUT_WORKDIR/Aligned.out.sam
mv "$CONT_OUTPUT_WORKDIR/Aligned.out.sam" "$CONT_OUTPUT_ALIGNED_FILE"
