#!/bin/bash

if (($# == 2)); then
    fragment_file=$1
    strand=$2
else
    echo "requires [fragment_file] [strand={-,+}]" 1>&2
    exit 1
fi

if [ "$strand" = "+" ]; then
    awk '{
        print $1,$2,$2+1,$4,NR,"+"
    }' $fragment_file | tr " " "\t"
elif [ "$strand" = "-" ]; then
    awk '{
        print $1,$3,$3+1,$4,NR,"-"
    }' $fragment_file | tr " " "\t"
else
    echo "stand must be either + or -" 1>&2
    exit 1
fi