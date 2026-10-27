#!/bin/bash

MYARR=$(awk '{print $2}' $1/metadata/CID_files.txt)
echo $MYARR
for i in ${MYARR[@]}
do
	if [ ! -f $1/RAW/$i ] && [ -f $1/RAW/$i.gz ]; then gunzip $1/RAW/$i.gz; fi
	if [ ! -f $1/RAW/$i ] && [ -f $1/RAW/$i.xz ]; then unxz -T16 $1/RAW/$i.xz; fi
	docker run --rm -v /public/:/public/ -w `pwd` mono mono ThermoRawFileParser/ThermoRawFileParser.exe -i=$1/RAW/$i -m=0 -o=$1/metadata/
done
