#!/bin/bash
# usage: bash get-rawfile-metadata.sh <PXD> [n_jobs (default: 4)]
NJOBS=${2:-4}

# on Ctrl+C or kill, stop the running jobs: docker passes the signal on to its container
trap 'kill $(jobs -p) 2>/dev/null; wait; exit 1' INT TERM

echo `pwd`
echo start $1
# write into a temporary folder and rename it to metadata only at the end,
# so an interrupted run does not leave an incomplete metadata folder behind
rm -rf "$1/metadata.tmp"
mkdir "$1/metadata.tmp"
# one ThermoRawFileParser container per RAW file, at most NJOBS at the same time
for f in $(find $1/RAW -maxdepth 1 -type f -iname '*.raw')
do
	docker run --rm -v /public/:/public/ -w `pwd` mono mono ThermoRawFileParser/ThermoRawFileParser.exe -i=$f -o=$1/metadata.tmp -m=0 &
	while [ $(jobs -rp | wc -l) -ge $NJOBS ]; do wait -n; done
done
wait
mv -T "$1/metadata.tmp" "$1/metadata"
echo end $1
