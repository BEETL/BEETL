#!/bin/bash

BEETL_INSTALL=/home/acox/Indel/tony/BEETL/

# Tidy up from prev runs

/bin/rm -rf bwt_phi bwt_phi_mut *-[ABC][0-9][0-9] phi_mut* phi_reads* mut_vs_* phi_rev.fa  genome_hits.txt

#exit

# Add mutations to phi genome
cat phi.fa | perl -alne '{ unless (defined($a)) { print; $a="" } else { $a.=$_ }} END { substr($a,4000,1)=""; substr($a,3000,1)="T"; substr($a,2000,0)="AA"; for ($i=0;$i<length($a);$i+=70) { print substr($a,$i,70) }}' > phi_mutated.fa

# Make read sets
cat phi.fa | tail -n+2 | perl -alne '{ $a.=$_ } END { for ($i=0;$i<length($a)-15;$i+=1) { print substr($a,$i,15) }}' > phi_reads_fwd.txt
cat phi_mutated.fa | tail -n+2 | perl -alne '{ $a.=$_ } END { for ($i=0;$i<length($a)-15;$i+=1) { print substr($a,$i,15) }}' > phi_mut_reads_fwd.txt
cat phi_reads_fwd.txt | perl -lane '{ $a=reverse($_); $a=~tr/ACGT/TGCA/; print $a }' > phi_reads_rev.txt
cat phi_mut_reads_fwd.txt | perl -lane '{ $a=reverse($_); $a=~tr/ACGT/TGCA/; print $a }' > phi_mut_reads_rev.txt

cat phi_reads_fwd.txt phi_reads_rev.txt > phi_reads.txt
cat phi_mut_reads_fwd.txt phi_mut_reads_rev.txt > phi_mut_reads.txt

# Make BWT of genome + reverse

cat phi.fa | perl -lane '{ if (defined($a)) { $a.=$_; } else { print; $a=""; } } END { $a=reverse($a); $a=~tr/ACGT/TGCA/; for ($i=0;$i<length($a);$i+=70) { print substr($a,$i,70) }}' > phi_rev.fa

${BEETL_INSTALL}/makeBWTSkew phi.fa bwt_phi.fa
${BEETL_INSTALL}/makeBWTSkew phi_rev.fa bwt_phi_rev.fa

${BEETL_INSTALL}/mergeBWT 0 phi_merge phi.fa phi_rev.fa
${BEETL_INSTALL}/mergeBWT 1 phi_merge phi.fa phi_rev.fa
${BEETL_INSTALL}/mergeBWT 2 phi_merge phi.fa phi_rev.fa
${BEETL_INSTALL}/mergeBWT 3 phi_merge phi.fa phi_rev.fa
${BEETL_INSTALL}/mergeBWT 4 phi_merge phi.fa phi_rev.fa
${BEETL_INSTALL}/mergeBWT 5 phi_merge phi.fa phi_rev.fa

# Make BWT of read sets

mkdir bwt_phi
cd bwt_phi
${BEETL_INSTALL}/BCRext ../phi_reads.txt
cd ..

mkdir bwt_phi_mut
cd bwt_phi_mut
${BEETL_INSTALL}/BCRext ../phi_mut_reads.txt
cd ..

# Compare read sets
${BEETL_INSTALL}/countWords 15 1 bwt_phi_mut/bwt-B?? bwt_phi/bwt-B0? > mut_vs_phi.txt

# Compare mutated read versus reference
${BEETL_INSTALL}/countWords -ref 15 1 bwt_phi_mut/bwt-B?? phi_merge-B?? > mut_vs_genome.txt

cat mut_vs_genome.txt  | ./find3.pl | sort -k 3,3n > genome_hits.txt
