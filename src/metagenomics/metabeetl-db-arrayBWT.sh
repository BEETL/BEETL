#!/bin/sh

#$ -cwd
#$ -S /bin/bash

SEQ_INPUT_DIR=/illumina/scratch/BWT/metagenomics/all.fna.tar/genomesToSingleSeq.out
BIN_DIR=/home/oschulz-trieglaff/code/BEETL_github/BEETL_DEV/scripts/metagenomics/database

echo "SGE_TASK_ID=${SGE_TASK_ID}"

cd ${SEQ_INPUT_DIR}
# build BWT of forward sequence
echo ${BIN_DIR}/metabeetl-db-makeBWTSkew ${SEQ_INPUT_DIR}/G_${SGE_TASK_ID} G_${SGE_TASK_ID}
${BIN_DIR}/metabeetl-db-makeBWTSkew ${SEQ_INPUT_DIR}/G_${SGE_TASK_ID} G_${SGE_TASK_ID}
# build BWT of reverse sequence
echo ${BIN_DIR}/metabeetl-db-makeBWTSkew ${SEQ_INPUT_DIR}/G_${SGE_TASK_ID}_rev G_${SGE_TASK_ID}_rev
${BIN_DIR}/metabeetl-db-makeBWTSkew ${SEQ_INPUT_DIR}/G_${SGE_TASK_ID}_rev G_${SGE_TASK_ID}_rev
