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

runCommand()
{
#COMMAND="${INSTALL}/Beetl ext -i ${TEST_FILE_SEQ} -p testBCRext -a -s" 
echo $0: Command to run: $COMMAND
#${TIME} $COMMAND 
${TIME} $COMMAND > $1.stdout 2>$1.stderr
echo $0: done
}

nextTest()
{
TEST_NUM=$(($TEST_NUM+1))
TEST_NAME=TEST_${TEST_NUM}
echo $0: $TEST_NAME - $1: `date`
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

TIME="/usr/bin/time"

if [ ! -e ${TIME} ]
then
	echo "$0: Can't find ${TIME}, so won't use it!"
	TIME=
else
	TIME="${TIME} -v"

fi

DIRS="BCR BCRext SAP"

READ_LENGTH_SHORT=37
TEST_FILE_SEQ=${PWD}/test.seq
TEST_FILE_FASTA=${PWD}/test.fasta
TEST_FILE_FASTQ=${PWD}/test.fastq
TEST_FILE_SEQ_SHORT=${PWD}/test_${READ_LENGTH_SHORT}.seq

echo $0: Cleaning up from previous test runs : `date`

/bin/rm -rf ${DIRS} ${TEST_FILE_SEQ} ${TEST_FILE_FASTA} ${TEST_FILE_FASTQ} ${TEST_FILE_SEQ_SHORT}

mkdir ${DIRS}

echo $0: Building test data files : `date`

# Definition of $c ensures 1 in 100 chars is an N

${PERL} -e '{ srand(1103); $c='ACGT' x 100 . 'NNNN' ; for ($i=0;$i<20000;$i++) { for ($j=0;$j<100;$j++) { print substr($c,rand(length($c)),1) } print "\n" }}' > ${TEST_FILE_SEQ}

cat ${TEST_FILE_SEQ} | ${PERL} -lane '{ $i++; print ">seq_${i}\n$_" }' > ${TEST_FILE_FASTA}

cat ${TEST_FILE_SEQ} | ${PERL} -lane '{ $i++; print "\@seq_${i}\n$_\n+"; print "X" x length($_) }' > ${TEST_FILE_FASTQ}

cat ${TEST_FILE_SEQ} | ${PERL} -lane "{ print substr(\$_,0, ${READ_LENGTH_SHORT}) }" > ${TEST_FILE_SEQ_SHORT}


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

TEST_NUM=0

nextTest "running BEETL BCRext code with ASCII output"

cd BCRext

COMMAND="${INSTALL}/Beetl ext -i ${TEST_FILE_SEQ} -p testBCRext -a -s" 
runCommand $TEST_NAME

cd ..

nextTest "Running BEETL BCRext code with run-length output"

cd BCRext

COMMAND="${INSTALL}/Beetl ext -i ${TEST_FILE_SEQ} -p testBCRextRunLength -a -s"
runCommand $TEST_NAME

cd ..

nextTest "Running BEETL BCRext code with Huffman output"

cd BCRext

COMMAND="${INSTALL}/Beetl ext -i ${TEST_FILE_SEQ} -p testBCRextHuffman -h -s"
runCommand $TEST_NAME

cd ..

nextTest "Partial comparison of BCRext output files"

cd BCRext

compareFiles testBCRext-B00 testBCRextRunLength-B00
compareFiles testBCRext-B00 testBCRextHuffman-B00

cd ..

nextTest "running BEETL BCRext code with FASTA input"

cd BCRext

COMMAND="${INSTALL}/Beetl ext -i ${TEST_FILE_FASTA} -p testBCRextFASTA -a" 
runCommand $TEST_NAME

cd ..

nextTest "Comparing BWT files for BCRext raw and FASTA input"

for d in 0 1 2 3 4 5
do 
  compareFiles BCRext/testBCRext-B0${d} BCRext/testBCRextFASTA-B0${d}
done



# Run BEETL BCR code

nextTest "Running BEETL BCR code on FASTA format file"

cd BCR

COMMAND="${INSTALL}/Beetl bcr -i ${TEST_FILE_FASTA} -o testBCR"
runCommand $TEST_NAME

cd ..

nextTest "Running BEETL BCR code on raw sequence file format"

cd BCR

COMMAND="${INSTALL}/Beetl bcr -i ${TEST_FILE_SEQ} -o testBCRSeq"
runCommand $TEST_NAME

cd ..

nextTest "Running BEETL BCR code on FASTQ file format"

cd BCR

COMMAND="${INSTALL}/Beetl bcr -i ${TEST_FILE_FASTQ} -o testBCRFASTQ"
runCommand $TEST_NAME

cd ..

nextTest "Comparing BWT files for BCR FASTA and raw sequence formats"

for d in 0 1 2 3 4 5
do 
  compareFiles BCR/testBCRSeq${d} BCR/testBCR${d}
done

nextTest "Comparing BWT files for BCR FASTA and FASTQ formats"

for d in 0 1 2 3 4 5
do 
  compareFiles BCR/testBCRFASTQ${d} BCR/testBCR${d}
done


nextTest "Comparing BWT files for BCRext and BCR"

for d in 0 1 2 3 4 5
do 
  compareFiles BCRext/testBCRext-B0${d} BCR/testBCR${d}
done



# Make BCR do BWT inversion 

nextTest "Running BEETL BCR inversion code"

cd BCR

#for d in testBCR?; do cp ${d} test_invert${d}; done
#touch test_invert
COMMAND="${INSTALL}/Beetl bcr -i testBCR -o testInvert.fa -m 1"
runCommand $TEST_NAME

nextTest "Checking BCR inversion against original output"

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

nextTest "Running BEETL BCRext code in SAP mode"

cd SAP

COMMAND="${INSTALL}/Beetl ext -i ${TEST_FILE_SEQ} -p testBCRext -a -sap -s"
runCommand $TEST_NAME

# Make BCR do BWT inversion 

nextTest "Running BEETL BCR inversion code on SAP-permuted BWT"

#cd SAP


#for d in BCR?; do cp ${d} test_invert${d}; done
touch testBCRext-B0
COMMAND="${INSTALL}/Beetl bcr -i testBCRext-B0 -o testInvert.fa -m 1"
runCommand $TEST_NAME

nextTest "Checking BCR inversion against original output"

grep -v '>' testInvert.fa | sort > testInvertSorted.txt

sort ${TEST_FILE_SEQ} > testFileSeqSorted.txt

compareFiles testFileSeqSorted.txt testInvertSorted.txt

cd ..

nextTest "running BEETL BCRext code on reads of length ${READ_LENGTH_SHORT}"

cd BCRext

COMMAND="${INSTALL}/Beetl ext -i ${TEST_FILE_SEQ_SHORT} -p testBCRext${READ_LENGTH_SHORT} -a -s" 
runCommand $TEST_NAME

cd ..

nextTest "Running BEETL BCR code on reads of length ${READ_LENGTH_SHORT}"

cd BCR

COMMAND="${INSTALL}/Beetl bcr -i ${TEST_FILE_SEQ_SHORT} -o testBCR${READ_LENGTH_SHORT}_"
runCommand $TEST_NAME

cd ..

nextTest "Comparing BWT files for BCRext and BCR on reads of length ${READ_LENGTH_SHORT}"

for d in 0 1 2 3 4 5
do 
  compareFiles BCRext/testBCRext${READ_LENGTH_SHORT}-B0${d} BCR/testBCR${READ_LENGTH_SHORT}_${d}
done

nextTest "Running BEETL BCR inversion code on reads of length ${READ_LENGTH_SHORT}"

cd BCR

#for d in testBCR?; do cp ${d} test_invert${d}; done
#touch test_invert
COMMAND="${INSTALL}/Beetl bcr -i testBCR${READ_LENGTH_SHORT}_ -o testInvert${READ_LENGTH_SHORT}.fa -m 1"
runCommand $TEST_NAME

nextTest "Checking BCR inversion against original output"

grep -v '>' testInvert${READ_LENGTH_SHORT}.fa > testInvert${READ_LENGTH_SHORT}.txt

compareFiles ${TEST_FILE_SEQ_SHORT} testInvert${READ_LENGTH_SHORT}.txt



# All tests done

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

