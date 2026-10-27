#!/bin/bash

# usage: bash Download-from-PRIDE-3.sh <ftp page of the project to download> <where to download>
# example: bash Download-from-PRIDE-3.sh https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/08/PXD010776/ /public/conode54_pride/PRIDE_DATA/PXD010776

cd $2
re="href=\"(\S+\.raw)\""
LINES=$(curl $1)
FILES=$(for l in ${LINES[@]}; do if [[ $l =~ $re ]]; then echo ${BASH_REMATCH[1]}; fi; done)
for f in ${FILES[@]}
do
    echo $f
done
for f in ${FILES[@]}
do
    wget $1/$f
done