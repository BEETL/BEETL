#!/bin/bash

## Copyright (c) 2011 Illumina, Inc.
##
## 
## This software is covered by the "Illumina Non-Commercial Use Software
## and Source Code License Agreement" and any user of this software or
## source file is bound by the terms therein (see accompanying file
## Illumina_Non-Commercial_Use_Software_and_Source_Code_License_Agreement.pdf)
##

## test_v_0_0_2.sh
## Unit tests for BEETL library code
## Author: Anthony J. Cox



compareFiles()
{
if [ ! -e "$1" ]
then 
    echo "$0: TEST FAILED! file $1 does not exist!"
    TEST_FAILED=1
fi

if [ ! -e "$2" ]
then 
    echo "$0: TEST FAILED! file $2 does not exist!"
    TEST_FAILED=1	
fi




a=`cat $1 | md5sum`
b=`cat $2 | md5sum`

#echo ${a} ${b}

if echo ${a} | grep -q ${b}
then 
    echo "$0: OK -  md5s for $1 and $2 coincide"
else 
    echo "$0: TEST FAILED! md5s for $1 and $2 differ!"
    TEST_FAILED=1
fi
}



PERL=/usr/bin/perl 

if [ ! -e ${PERL} ]
then
	echo "$0: A Perl installation is needed to run this script, can't find one at ${PERL}."
	echo "$0: Either install Perl, or (more likely) modify variable PERL in $0 to point to a Perl installation on your system."
	exit -1
fi


INSTALL=`dirname $0`

if [ $INSTALL = "." ]
then
	echo $0: "Running this script will create new subdirectories of the current working directory for test purposes."
	echo $0: "Your current working directory seems to be the BEETL installation directory."
	echo $0: "To test BEETL, please move to a new directory (e.g. /tmp) and re-run this script."
	exit -1
fi

if [ $INSTALL = $PWD ]
then
	echo $0: "Running this script will create new subdirectories of the current working directory for test purposes."
	echo $0: "Your current working directory seems to be the BEETL installation directory."
	echo $0: "To test BEETL, please move to a new directory (e.g. /tmp) and re-run this script."
	exit -1
fi

TIME="/usr/bin/time -v"

DIRS="v1_30 BCR BCRext SAP"


TEST_FILE_SEQ=${PWD}/test.seq
TEST_FILE_FASTA=${PWD}/test.fasta
TEST_FILE_FASTQ=${PWD}/test.fastq

echo $0: Cleaning up from previous test runs : `date`

/bin/rm -rf ${DIRS} ${TEST_FILE_SEQ} ${TEST_FILE_FASTA} ${TEST_FILE_FASTQ}

mkdir ${DIRS}

echo $0: Building test data files : `date`

# Definition of $c ensures 1 in 100 chars is an N

${PERL} -e '{ srand(1103); $c='ACGT' x 100 . 'NNNN' ; for ($i=0;$i<20000;$i++) { for ($j=0;$j<100;$j++) { print substr($c,rand(length($c)),1) } print "\n" }}' > ${TEST_FILE_SEQ}

cat ${TEST_FILE_SEQ} | ${PERL} -lane '{ $i++; print ">seq_${i}\n$_" }' > ${TEST_FILE_FASTA}

cat ${TEST_FILE_SEQ} | ${PERL} -lane '{ $i++; print "\@seq_${i}\n$_\n+"; print "X" x length($_) }' > ${TEST_FILE_FASTQ}

#exit;

# Run CPM 2011 code

# Don't do this - not portable - if generated BWTs invert to original sequences then that is a good enough test

#echo $0: Running CPM version of BCRext code : `date`
#
#cd v1_30
#
#${TIME} /home/acox/Indel/tony/Misc/embwt_v1_30 ${TEST_FILE_SEQ} >  run.stdout 2>run.stderr
#
#cd ..

# Run BEETL BCRext code

echo $0: Running BEETL BCRext code with ASCII output: `date`

cd BCRext

${TIME} ${INSTALL}/Beetl ext -i ${TEST_FILE_SEQ} -p testBCRext -a  >  run.stdout 2>run.stderr

cd ..

echo $0: Running BEETL BCRext code with run-length output: `date`

cd BCRext

${TIME} ${INSTALL}/Beetl ext -i ${TEST_FILE_SEQ} -p testBCRextRunLength -a  >  run_rle.stdout 2>run_rle.stderr

cd ..

echo $0: Running BEETL BCRext code with Huffman output: `date`

cd BCRext

${TIME} ${INSTALL}/Beetl ext -i ${TEST_FILE_SEQ} -p testBCRextHuffman -h  >  run_huff.stdout 2>run_huff.stderr

cd ..

echo $0: Partial comparison of BCRext output files: `date`

cd BCRext

compareFiles testBCRext-B00 testBCRextRunLength-B00
compareFiles testBCRext-B00 testBCRextHuffman-B00

cd ..

# Run BEETL BCR code

echo $0: Running BEETL BCR code on FASTA format file: `date`

cd BCR

${TIME} ${INSTALL}/Beetl bcr -i ${TEST_FILE_FASTA} -o testBCR  >  run.stdout 2>run.stderr

cd ..

echo $0: Running BEETL BCR code on raw sequence file format: `date`

cd BCR

${TIME} ${INSTALL}/Beetl bcr -i ${TEST_FILE_SEQ} -o testBCRSeq  >  run.stdout 2>run.stderr

cd ..

echo $0: Comparing BWT files for BCR FASTA and raw sequence formats: `date`

for d in 0 1 2 3 4 5
do 
  compareFiles BCR/testBCRSeq${d} BCR/testBCR${d}
done
cd ..

echo $0: Comparing BWT files for BCRext and BCR: `date`

for d in 0 1 2 3 4 5
do 
  compareFiles BCRext/testBCRext-B0${d} BCR/testBCR${d}
done




# Make BCR do BWT inversion 

echo $0: Running BEETL BCR inversion code : `date`

cd BCR

#for d in testBCR?; do cp ${d} test_invert${d}; done
#touch test_invert
${TIME} ${INSTALL}/Beetl bcr -i testBCR -o testInvert.fa -m 1 > invert.stdout 2> invert.stderr

echo $0: Checking BCR inversion against original output : `date`


#for d in cyc.*.txt
#do
#  CYCLE=`echo ${d} | sed 's/^cyc.//' | sed 's/.txt$//'`
##  echo ${CYCLE}
#  cat ${TEST_FILE_SEQ} | perl -ane "{ print substr(\$_,${CYCLE},1) }" > check_${CYCLE}.txt
#  compareFiles cyc.${CYCLE}.txt check_${CYCLE}.txt
#done
# BCR now does the restriping into a fasta file so this bit no longer needed

grep -v '>' testInvert.fa > testInvert.txt

compareFiles ${TEST_FILE_SEQ} testInvert.txt

cd ..

# Run BEETL BCRext code in SAP mode

echo $0: Running BEETL BCRext code in SAP mode: `date`

cd SAP

${TIME} ${INSTALL}/Beetl ext -i ${TEST_FILE_SEQ} -p testBCRext -a -sap  >  run.stdout 2>run.stderr

# Make BCR do BWT inversion 

echo $0: Running BEETL BCR inversion code on SAP-permuted BWT: `date`

#cd SAP


#for d in BCR?; do cp ${d} test_invert${d}; done
touch testBCRext-B0
${TIME} ${INSTALL}/Beetl bcr -i testBCRext-B0 -o testInvert.fa -m 1 > invert.stdout 2> invert.stderr

echo $0: Checking BCR inversion against original output : `date`


grep -v '>' testInvert.fa | sort > testInvertSorted.txt

sort ${TEST_FILE_SEQ} > testFileSeqSorted.txt

compareFiles testFileSeqSorted.txt testInvertSorted.txt


if [ -n "$TEST_FAILED" ]
then 
    echo "$0: ***                         *** "
    echo "$0: *** Test failures reported! *** "
    echo "$0: ***                         *** "
    echo $0: Test run finished : `date`
    exit -1;
else
    echo "$0: ***                           *** "
    echo "$0: *** No test failures reported *** "
    echo "$0: ***                           *** "
    echo $0: Test run finished : `date`
    exit 0;
fi

