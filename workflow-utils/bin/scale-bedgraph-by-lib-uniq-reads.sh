#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <bedgraph file> <read count file>"
    exit 1
fi

# compute library size of fragment-extended bedgraph
lib_size=$(cat $2)
scale_factor=$(bc <<< "scale=4; 100000000/${lib_size}")

# scale bedgraph
awk -v scale_factor=$scale_factor 'BEGIN {FS="\t"; OFS="\t"}; {print $1, $2, $3, $4 * scale_factor}' $1