#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <bowtie log>"
    exit 1
fi

reads_processed=$(grep 'reads processed' $1 | cut -f 4 -d ' ')
reads_mapped=$(grep 'Reported' $1 | cut -f 2 -d ' ')
echo -e $reads_processed"\t"$reads_mapped
