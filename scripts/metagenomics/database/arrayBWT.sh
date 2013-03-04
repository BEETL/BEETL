#!/bin/sh

$-cwd

cd /illumina/scratch/BWT/metagenomics/allNucleotides/workDir/singleSeq
/home/cander/Beetl/mergeChromosoms/makeBWTSkew /illumina/scratch/BWT/metagenomics/allNucleotides/workDir/singleSeq/G_"$SGE_TASK_ID"_rev G_"$SGE_TASK_ID"_rev

echo /home/cander/Beetl/mergeChromosoms/makeBWTSkew /illumina/scratch/BWT/metagenomics/allNucleotides/workDir/singleSeq/G_'$SGE_TASK_ID'_rev G_'$SGE_TASK_ID'_rev