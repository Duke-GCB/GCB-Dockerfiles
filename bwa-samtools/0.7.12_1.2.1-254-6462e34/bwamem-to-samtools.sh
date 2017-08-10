#!/bin/bash

# Shell script to run bwa mem and pipe output to samtools
# luckily, bwa outputs to stdout, so this tool can too

bwa mem $@ | samtools view -bS -

