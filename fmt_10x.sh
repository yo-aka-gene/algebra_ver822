#!/bin/bash

CMDNAME=`basename $0`
if [ $# -ne 1 ]; then
    echo "Error: $CMDNAME expects single argument; e.g., $CMDNAME dirname"
    exit 1
fi

chmod 777 $1/*.gz
for i in $(ls $1/*_barcodes.tsv.gz); do mkdir $1/$(basename $i _barcodes.tsv.gz); done
for i in $(ls $1/*_barcodes.tsv.gz); do mv $i $1/$(basename $i _barcodes.tsv.gz)/barcodes.tsv.gz; done
for i in $(ls $1/*_features.tsv.gz); do mv $i $1/$(basename $i _features.tsv.gz)/features.tsv.gz; done
for i in $(ls $1/*_matrix.mtx.gz); do mv $i $1/$(basename $i _matrix.mtx.gz)/matrix.mtx.gz; done

exit 0
