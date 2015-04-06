#!/bin/bash

set -e

# This is a script for a docker image to run trimmomatic, a java program
# See http://www.usadellab.org/cms/?page=trimmomatic
#
# CONT prefix means this variable is expected to be set in the container at runtime

# INPUTS
#
# CONT_INPUT_FASTQ_READ_1 - Required for both Paired end and single end
# CONT_INPUT_FASTQ_READ_2 - Required only for Paired end
# CONT_INPUT_ADAPTERS_FILE - Required
#
# PARAMETERS
#
# CONT_PARAM_TRIM_MODE  PE or SE, required.
# CONT_PARAM_QUALITY_SCORE  default -phred33
# CONT_PARAM_ILLUMINACLIP_SETTINGS - settings for <seed mismatches>:<palindrome clip threshold>:<simple clip threshold>
# CONT_PARAM_ADDITIONAL_STEPS default LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30
#
# OUTPUTS
#
# CONT_OUTPUT_SE_TRIMMED - Required for SE
# CONT_OUTPUT_PE_TRIMMED_P1 - Required for PE
# CONT_OUTPUT_PE_TRIMMED_P2 - Required for PE
# CONT_OUTPUT_PE_TRIMMED_U1 - Required for PE
# CONT_OUTPUT_PE_TRIMMED_U2 - Required for PE

# Check that trim mode is set
[ -z "$CONT_PARAM_TRIM_MODE" ] && echo "Error: The CONT_PARAM_TRIM_MODE variable must be set" && exit 1

# Check parameters based on trim mode
if [ "PE" -eq "$CONT_PARAM_TRIM_MODE" ]
then
  # If PE, need 2 inputs and 4 outputs
  for VARNAME in CONT_INPUT_FASTQ_READ_1 CONT_INPUT_FASTQ_READ_2 CONT_INPUT_ADAPTERS_FILE do
    eval INFILE=\$$VARNAME
    # Check that input variables are set
    [ -z "$INFILE" ] && echo "Error: The $VARNAME variable must be set" && exit 1
    # Check that input files are readable
    [ ! -r "$INFILE" ] && echo "Error: input file $INFILE does not exist or cannot be read" && exit 1
  done
  for VARNAME in CONT_OUTPUT_PE_TRIMMED_P1 CONT_OUTPUT_PE_TRIMMED_P2 CONT_OUTPUT_PE_TRIMMED_U1 CONT_OUTPUT_PE_TRIMMED_U2 do
    eval OUTFILE=\$$VARNAME
    # Check that output variables are set
    [ -z "$OUTFILE" ] && echo "Error: The $VARNAME variable must be set" && exit 1
    # Check that output directories are writable
    [ ! -w $(dirname "$OUTFILE" ) ] && echo "Error: output file $OUTFILE is not writable" && exit 1
  done
  TRIMMOMATIC_BIN=`which TrimmomaticPE`
  $INPUTS="$CONT_INPUT_FASTQ_READ_1 $CONT_INPUT_FASTQ_READ_2"
  $OUTPUTS="$CONT_OUTPUT_PE_TRIMMED_P1 $CONT_OUTPUT_PE_TRIMMED_P2 $CONT_OUTPUT_PE_TRIMMED_U1 $CONT_OUTPUT_PE_TRIMMED_U2"
elif [ "SE" -eq "$CONT_PARAM_TRIM_MODE" ]
then
  # If SE, one input and 1 output
  for VARNAME in CONT_INPUT_FASTQ_READ_1 CONT_INPUT_ADAPTERS_FILE do
    eval INFILE=\$$VARNAME
    # Check that input variables are set
    [ -z "$INFILE" ] && echo "Error: The $VARNAME variable must be set" && exit 1
    # Check that input files are readable
    [ ! -r "$INFILE" ] && echo "Error: input file $INFILE does not exist or cannot be read" && exit 1
  done
  # Check that output variables are set
  [ -z "$CONT_OUTPUT_SE_TRIMMED" ] && echo "Error: The CONT_OUTPUT_SE_TRIMMED variable must be set" && exit 1
  # Check that output directories are writable
  [ ! -w $(dirname "$CONT_OUTPUT_SE_TRIMMED" ) ] && echo "Error: output file $CONT_OUTPUT_SE_TRIMMED is not writable" && exit 1
  TRIMMOMATIC_BIN=`which TrimmomaticSE`
  $INPUTS=$CONT_INPUT_FASTQ_READ_1
  $OUTPUTS=$CONT_OUTPUT_SE_TRIMMED
else
  echo "Error: Invalid CONT_PARAM_TRIM_MODE: $CONT_PARAM_TRIM_MODE" && exit 1
fi

# Mode, input and output set. Check other parameters
# If non-default values provided to container, use them.
# ILLUMINACLIP_STEP includes the adapters file
ILLUMINACLIP_STEP="ILLUMINACLIP:$CONT_INPUT_ADAPTERS_FILE:2:30:15"
ADDITIONAL_STEPS="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30"
THREADS="-threads 1"
QUALITY_SCORE="-phred33"

if [ ! -z "$CONT_PARAM_ILLUMINACLIP_SETTINGS" ]; then
  ILLUMINACLIP_STEP="ILLUMINACLIP:$CONT_INPUT_ADAPTERS_FILE:$CONT_PARAM_ILLUMINACLIP_SETTINGS"
fi

if [ ! -z "$CONT_PARAM_ADDITIONAL_STEPS" ]; then
  ADDITIONAL_STEPS=$CONT_PARAM_ADDITIONAL_STEPS
fi

if [ ! -z "$CONT_PARAM_QUALITY_SCORE" ]; then
  $QUALITY_SCORE=$CONT_PARAM_QUALITY_SCORE
fi

TRIMMOMATIC_CMD="$TRIMMOMATIC_BIN $THREADS $QUALITY_SCORE $INPUTS $OUTPUTS $ILLUMINACLIP_STEP $ADDITIONAL_STEPS"
TRIMMOMATIC_VERSION=$(readlink "/usr/share/java/trimmomatic.jar")

echo
echo "Starting $0..."
echo "$TRIMMOMATIC_VERSION"
echo "$TRIMMOMATIC_CMD"
sh -c "$TRIMMOMATIC_CMD"
