#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Requires [FASTA] [BIAS]" 1>&2
    exit 1
fi

faidx $1 --bed $2 --transform nucleotide --lazy -s Q | cut -f 4-8 | \
    sed '1d' | awk '{print ($4 + $3)/($1+$2+$3+$4), $1+$2+$3+$4+$5}' | tr " " "\t"
