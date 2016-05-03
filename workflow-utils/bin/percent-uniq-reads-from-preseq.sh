#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <preseq c_curve outfile>"
    exit 1
fi

total='15000000'
distinct=$(grep $total $1 | cut -f 2 | awk '{ print +$1 }' )
if [[ $distinct =~ ^\ *$ ]]; then
    perc_distinct="NA"
else
    perc_distinct=$(bc <<< "scale=4; $distinct/$total")
fi
echo $perc_distinct