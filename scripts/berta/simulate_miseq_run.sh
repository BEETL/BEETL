#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3"  ] ; then
	echo "Usage: $0 <INPUT_DIR> <OUTPUT_DIR> <NUM_CYCLES>";
  exit 1
fi

if ! [ -d "$1" ]; then
	echo "Directory $1 does not seem to exist.";
  exit 1
fi  

if ! [ -d "$2" ]; then
	echo "Directory $2 does not seem to exist.";
  exit 1
fi  

INDIR=$1;
OUTDIR=$2;
NUM_CYCLES=$3;

#MISEQ_CYCLE_TIME=2s # for testing
MISEQ_CYCLE_TIME=6m

BCL_SUFFIX="Data/Intensities/BaseCalls/L001"

BCL_IN_DIR="$INDIR/$BCL_SUFFIX"
BCL_OUT_DIR="$OUTDIR/$BCL_SUFFIX"

echo "Reading from $BCL_IN_DIR";
echo "Writing to $BCL_OUT_DIR";
echo "Simulating $NUM_CYCLES cycles";

if ! [ -a "$BCL_OUT_DIR" ]; then
	echo "$BCL_OUT_DIR does not seem to exist. Creating it.";
	mkdir -p $BCL_OUT_DIR;
fi

for i in `seq 1 $NUM_CYCLES`
do
	echo "Simulating cycle $i";
	CURRENT_BCL="${BCL_IN_DIR}/C${i}.1"
	echo ln -s ${CURRENT_BCL} $BCL_OUT_DIR
	ln -s ${CURRENT_BCL} $BCL_OUT_DIR
	echo "Sleeping $MISEQ_CYCLE_TIME";
	sleep $MISEQ_CYCLE_TIME;
done

