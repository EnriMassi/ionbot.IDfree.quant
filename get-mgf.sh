#!/bin/bash
# usage: bash get-mgf.sh <PXD> [n_jobs (default: 4)]
NJOBS=${2:-4}

MYARR=$(awk -F '\t' '{print $2}' $1/metadata/HCD_files.txt)

# on Ctrl+C or kill, stop the running jobs: docker passes the signal on to its container
trap 'kill $(jobs -p) 2>/dev/null; wait; exit 1' INT TERM

# write into a temporary folder and rename it to MGF only at the end,
# so an interrupted run does not leave an incomplete MGF folder behind
rm -rf "$1/MGF.tmp"
mkdir "$1/MGF.tmp"
# decompress and convert each RAW file in a subshell, at most NJOBS at the same time
for i in ${MYARR[@]}
do
	(
	if [ ! -f $1/RAW/$i ] && [ -f $1/RAW/$i.gz ]; then gunzip $1/RAW/$i.gz; fi
	if [ ! -f $1/RAW/$i ] && [ -f $1/RAW/$i.xz ]; then unxz -T16 $1/RAW/$i.xz; fi
	exec docker run --rm -v /public/:/public/ -w `pwd` mono mono ThermoRawFileParser/ThermoRawFileParser.exe -i=$1/RAW/$i -o=$1/MGF.tmp -f=0
	) &
	while [ $(jobs -rp | wc -l) -ge $NJOBS ]; do wait -n; done
done
wait
mv -T "$1/MGF.tmp" "$1/MGF"
