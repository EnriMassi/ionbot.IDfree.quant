#!/bin/bash
echo start $1
MYARR=$(awk -F '\t' '{print $2}' $1/metadata/HCD_files.txt)
for i in ${MYARR[@]}
do
	if [ ! -f $1/RAW/$i ] && [ -f $1/RAW/$i.gz ]; then gunzip $1/RAW/$i.gz; fi
	if [ ! -f $1/RAW/$i ] && [ -f $1/RAW/$i.xz ]; then unxz -T16 $1/RAW/$i.xz; fi
done
echo end $1
