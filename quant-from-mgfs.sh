#!/bin/bash
# usage: bash quant-from-mgfs.sh [n_pxds (default: 2)] [n_jobs (default: 4)]
# n_pxds projects run at the same time, each with up to n_jobs RAW files in parallel,
# so up to n_pxds * n_jobs ThermoRawFileParser containers can run at once.
# FlashLFQ runs on one mzML at a time within each project.
NPXD=${1:-2}
NJOBS=${2:-4}

quantify() {
	echo "_______"
	echo $1
	python -u RAW-files-preproc.py -j $2 $1
	python -u ionbot.quant/mgf2FlashLFQ.py $1
	gzip -f $1/MGF/*.mgf
	gzip -f $1/mzml/*.mzML
	rm $1/FlashLfqSettings.toml $1/flashlfq.tsv  $1/QuantifiedPeaks.tsv
	echo "_______"
}
export -f quantify

# skip duplicate PXDs: two processes must never work on the same folder
MYARR=$(awk '!seen[$1]++ {print $1}' pxds.to.quantify.txt)
echo $MYARR
# the output of each project goes to <PXD>/quant.log
printf '%s\n' ${MYARR[@]} | xargs -P $NPXD -I{} bash -c 'echo start $1; quantify $1 $2 > $1/quant.log 2>&1; echo end $1' _ {} $NJOBS
