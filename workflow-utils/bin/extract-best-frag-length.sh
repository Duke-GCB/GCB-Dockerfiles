#!/bin/bash

frag=0;
for ii in $(cut -f3 $1 | tr "," " ");
do
    if [[ $ii -gt 0 ]]; then
        frag=$ii;
        break;
    fi;
done;
echo $frag;