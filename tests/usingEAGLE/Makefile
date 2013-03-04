TOP                 = $(abspath $(dir $(lastword $(MAKEFILE_LIST))))

EAGLE_PATH         ?= /illumina/development/EAGLE/EAGLE-1.4.0
BCL2FASTQ_PATH     ?= /illumina/development/casava/bcl2fastq-1.8.3b2
ISAAC_PATH         ?= /illumina/development/iSAAC/testing
STARLING2VCF_PATH  ?= /home/ljanin/workspace/WhiteFly/project/FeatureCallingEval/bin
BEETL_PATH          = $(TOP)/..
BEETL_LCP_PATH      = /home/ljanin/workspace/git/BEETL_DEV_for_LCP/BEETL_DEV


#TEST               = EColi
TEST               ?= HumanDiploidChr21
MAX_RAM            ?= 60

CUTOFF_THRESHOLD   ?= 19
CUTOFF_QUALITY     ?= 40
CUTOFF_TMP_QUAL     = 255


ifeq ($(TEST),EColi)
REFERENCE_DIR       = /illumina/development/iSAAC/iGenomes/Escherichia_coli_K_12_MG1655/NCBI/2001-10-15/Sequence/iSAACIndex.20120203
REFERENCE_GENOME    = ${REFERENCE_DIR}/GATK/genome.fa
ISAAC_REFERENCE     = ${REFERENCE_DIR}/genome.fa-32mer-2bit-SortedReference.xml
CHROM               = chr
COVERAGE_DEPTH     ?= 5
VARIANT_LIST	    = 
TEST_ID             = EColi_${COVERAGE_DEPTH}x

else ifeq ($(TEST),HumanDiploidChr21)
REFERENCE_GENOME    = /illumina/scratch/eagle/References/iGenomes_hg19_with_fai/chr21.fa
ISAAC_REFERENCE     = /illumina/scratch/eagle/References/forISAAC/HumanChr21/chr21.fa-32mer-6bit-SortedReference.xml
CHROM               = chr21
COVERAGE_DEPTH     ?= 15
VARIANT_LIST	    = /illumina/scratch/eagle/RealisticInputs/dbSnp_set2_UCSC_chr21.vcf
EXTRA_EAGLE_OPTIONS = --genome-mutator-options="--organism-ploidy=2"
TEST_ID             = HumanDiploidChr21_2x${COVERAGE_DEPTH}x

else
$(error Invalid test name)

endif


ifneq (,$(VARIANT_LIST))
VARIANT_LIST_OPTION = --variant-list=$(VARIANT_LIST)
endif
CUTOFF_QUALITY_FASTQ_OCTAL = $(shell printf "%03o" $$[${CUTOFF_QUALITY}+33])
CUTOFF_QUALITY_FASTQ_CHAR   = $(shell printf "\\${CUTOFF_QUALITY_FASTQ_OCTAL}" )
CUTOFF_TMP_QUAL_OCTAL      = $(shell printf "%03o" ${CUTOFF_TMP_QUAL})
CUTOFF_TMP_QUAL_CHAR       = $(shell printf "\\${CUTOFF_TMP_QUAL_OCTAL}" )

CUTOFF_LEVEL1_DIR = cutoff_threshold${CUTOFF_THRESHOLD}
CUTOFF_LEVEL2_DIR = ${CUTOFF_LEVEL1_DIR}/Q${CUTOFF_QUALITY}


.PHONY: all sim-pre sim fastq-pre fastq fasta bwt lcp cutoff unbwt pre-isaac post-isaac hist report
all: report20
sim-pre: EAGLE_${TEST_ID}/Makefile
sim: EAGLE_${TEST_ID}/.L001_1.bcl.completed
fastq-pre: FastQ/Makefile
fastq: reads.fastq
fasta: reads.fa
bwt: BWT/reads.bwt_bcr
lcp: LCP/reads.bwt_bcr.part
cutoff: ${CUTOFF_LEVEL1_DIR}/cutoffBWT/reads.bwt_bcr.part1.qual
unbwt: ${CUTOFF_LEVEL1_DIR}/decodedBWT.fastq
pre-isaac: pre_iSAAC/Aligned/Projects/default/default/sorted.bam
post-isaac: ${CUTOFF_LEVEL2_DIR}/post_iSAAC/Aligned/Projects/default/default/sorted.bam

report: ${CUTOFF_LEVEL2_DIR}/report.txt
	cat $^

report20: ${CUTOFF_LEVEL2_DIR}/report_aboveQ20.txt
	cat $^

forceBWT:
	rm -rf BWT
	$(MAKE) -f ${TOP}/Makefile bwt


RunInfo.xml: ${EAGLE_PATH}/share/EAGLE/RunInfo/RunInfo_PairedReads1x1Tiles.xml
	cp $^ $@
#	sed 's/101/100/' $^ > $@

#QualityTable.qval: /illumina/development/EAGLE/EAGLE-1.4.0/share/EAGLE-1.4.0/QualityTables/DefaultQualityTable.read1.length101.qval
#	grep -vP '^101' $^ > $@

EAGLE_${TEST_ID}/Makefile: RunInfo.xml #QualityTable.qval
	${EAGLE_PATH}/bin/configureEAGLE.pl \
	  --run-info=RunInfo.xml \
	  --reference-genome=${REFERENCE_GENOME} \
	  --coverage-depth=${COVERAGE_DEPTH} \
	  ${VARIANT_LIST_OPTION} \
	  ${EXTRA_EAGLE_OPTIONS} \
	  --run-folder=RunFolder \
	  EAGLE_${TEST_ID}
#	  --quality-table=QualityTable.qval \

EAGLE_${TEST_ID}/.L001_1.bcl.completed: EAGLE_${TEST_ID}/Makefile
	time $(MAKE) -C EAGLE_${TEST_ID}

FastQ/Makefile: EAGLE_${TEST_ID}/.L001_1.bcl.completed
	${BCL2FASTQ_PATH}/bin/configureBclToFastq.pl \
	  --use-bases-mask "y100n1,y100n1" \
	  --input-dir EAGLE_${TEST_ID}/RunFolder/Data/Intensities/BaseCalls \
	  --output-dir FastQ \
	; echo "FastQ configured"

FastQ/Project_FC1234TST/Sample_lane1/lane1_NoIndex_L001_R1_001.fastq.gz: FastQ/Makefile
	time $(MAKE) -C FastQ l1.done

reads.fastq: FastQ/Project_FC1234TST/Sample_lane1/lane1_NoIndex_L001_R1_001.fastq.gz
	zcat FastQ/Project_FC1234TST/Sample_lane*/lane*_NoIndex_L00*_R1_*.fastq.gz FastQ/Project_FC1234TST/Sample_lane*/lane*_NoIndex_L00*_R2_*.fastq.gz > $@

reads.fa: reads.fastq
	time sed -n 's/^/>/;p;n;p;n;n' $^ > $@

BWT/reads.bwt_bcr: reads.fastq
	mkdir -p BWT \
	&& cd BWT \
	&& time ${BEETL_PATH}/Beetl bcr -i ../$^ -o ../$@.part \
	&& cd .. \
	&& cat $@.part* > $@

LCP/reads.bwt_bcr.part: reads.fa
	mkdir -p LCP \
	&& cd LCP \
	&& time ${BEETL_LCP_PATH}/Beetl bcr -i ../$^ -o ../$@


define RSCRIPT
	for (i in 1:6) {\n
		png( paste( "hist_LCP", i, ".png", sep="" ) )\n
		f = file( paste( "LCP/reads.bwt_bcr.part", i, ".lcp", sep="" ), "rb")\n
		x <- readBin( f, integer(), n=100000000, size=4, endian="little" )\n
		hist( x, breaks=seq(min(x),max(x)) )\n
		dev.off()\n
	}
endef
export RSCRIPT

hist_LCP1.png: LCP/reads.bwt_bcr.part
	echo -e ""$$RSCRIPT"" > script.R \
	; R CMD BATCH script.R \
	; echo "hist generated"

hist: hist_LCP1.png
	gqview $^

#BWT/reads.bwt_ext: reads.fa
#	mkdir -p BWT \
#	&& time ../Beetl ext -i $^ -p $@.part -a \
#	&& cat $@.part-B* > $@

#decodedBWT.fa: BWT/reads.bwt_bcr
#	${BEETL_PATH}/Beetl bcr -m 1 -i $^.part -o $@

${CUTOFF_LEVEL1_DIR}/cutoffBWT/.sentinel:
	mkdir -p $(dir $@) \
	&& touch $@

${CUTOFF_LEVEL2_DIR}/.sentinel:
	mkdir -p $(dir $@) \
	&& touch $@

${CUTOFF_LEVEL1_DIR}/cutoffBWT/reads.bwt_bcr.part%.qual: BWT/reads.bwt_bcr LCP/reads.bwt_bcr.part ${CUTOFF_LEVEL1_DIR}/cutoffBWT/.sentinel
	ln -s ../../BWT/reads.bwt_bcr.part$* $(dir $@) \
	&& time ${BEETL_PATH}/scripts/lcp/applyLcpCutoff.pl -i BWT/reads.bwt_bcr.part$*.qual -l LCP/reads.bwt_bcr.part$*.lcp -b BWT/reads.bwt_bcr.part$* -o $@ -c ${CUTOFF_THRESHOLD} -r ${CUTOFF_TMP_QUAL} > $@.stdout

${CUTOFF_LEVEL1_DIR}/decodedBWT.fastq: ${CUTOFF_LEVEL1_DIR}/cutoffBWT/.sentinel ${CUTOFF_LEVEL1_DIR}/cutoffBWT/reads.bwt_bcr.part0.qual ${CUTOFF_LEVEL1_DIR}/cutoffBWT/reads.bwt_bcr.part1.qual ${CUTOFF_LEVEL1_DIR}/cutoffBWT/reads.bwt_bcr.part2.qual ${CUTOFF_LEVEL1_DIR}/cutoffBWT/reads.bwt_bcr.part3.qual ${CUTOFF_LEVEL1_DIR}/cutoffBWT/reads.bwt_bcr.part4.qual ${CUTOFF_LEVEL1_DIR}/cutoffBWT/reads.bwt_bcr.part5.qual ${CUTOFF_LEVEL1_DIR}/cutoffBWT/reads.bwt_bcr.part6.qual
	cd ${CUTOFF_LEVEL1_DIR}/cutoffBWT \
	&& time ${BEETL_PATH}/Beetl bcr -m 1 -i reads.bwt_bcr.part -o ../../$@

${CUTOFF_LEVEL2_DIR}/decodedBWT.fastq: ${CUTOFF_LEVEL1_DIR}/decodedBWT.fastq ${CUTOFF_LEVEL2_DIR}/.sentinel
	time cat $< | tr "\${CUTOFF_TMP_QUAL_OCTAL}" "${CUTOFF_QUALITY_FASTQ_CHAR}" > $@


#decodedBWT.fastq.gz: ${CUTOFF_LEVEL2_DIR}/decodedBWT.fastq
#	gzip -c $^ > $@

pre_iSAAC/inputs/.sentinel:
	mkdir -p $(dir $@) \
	&& touch $@

# % is for "reads"(.fastq) and "${CUTOFF_LEVEL2_DIR}/decodedBWT"(.fastq)
%.fastq.0: %.fastq
	let "LINECOUNT = $(shell wc -l $^ | cut -f 1 -d ' ' ) / 2" \
	; time split -d -a 1 -l $${LINECOUNT} $^ $^.

pre_iSAAC/inputs/lane1_read1.fastq.gz: pre_iSAAC/inputs/.sentinel reads.fastq.0
	gzip -c reads.fastq.0 > $@

pre_iSAAC/inputs/lane1_read2.fastq.gz: pre_iSAAC/inputs/.sentinel reads.fastq.0
	gzip -c reads.fastq.1 > $@

pre_iSAAC/Aligned/Projects/default/default/sorted.bam: pre_iSAAC/inputs/lane1_read1.fastq.gz pre_iSAAC/inputs/lane1_read2.fastq.gz
	time ${ISAAC_PATH}/bin/isaac-align \
	  -r ${ISAAC_REFERENCE} \
	  --base-calls-format fastq-gz \
	  --use-bases-mask="y100,y100" \
	  -b pre_iSAAC/inputs \
	  -o pre_iSAAC/Aligned \
	  -t pre_iSAAC/Temp \
	  -m ${MAX_RAM} &> pre_iSAAC/isaac.log

${CUTOFF_LEVEL2_DIR}/post_iSAAC/inputs/.sentinel:
	mkdir -p $(dir $@) \
	&& touch $@

${CUTOFF_LEVEL2_DIR}/post_iSAAC/inputs/lane1_read1.fastq.gz: ${CUTOFF_LEVEL2_DIR}/decodedBWT.fastq.0 ${CUTOFF_LEVEL2_DIR}/post_iSAAC/inputs/.sentinel
	gzip -c ${CUTOFF_LEVEL2_DIR}/decodedBWT.fastq.0 > $@

${CUTOFF_LEVEL2_DIR}/post_iSAAC/inputs/lane1_read2.fastq.gz: ${CUTOFF_LEVEL2_DIR}/decodedBWT.fastq.0 ${CUTOFF_LEVEL2_DIR}/post_iSAAC/inputs/.sentinel
	gzip -c ${CUTOFF_LEVEL2_DIR}/decodedBWT.fastq.1 > $@

${CUTOFF_LEVEL2_DIR}/post_iSAAC/Aligned/Projects/default/default/sorted.bam: ${CUTOFF_LEVEL2_DIR}/post_iSAAC/inputs/lane1_read1.fastq.gz ${CUTOFF_LEVEL2_DIR}/post_iSAAC/inputs/lane1_read2.fastq.gz
	time ${ISAAC_PATH}/bin/isaac-align \
	  -r ${ISAAC_REFERENCE} \
	  --base-calls-format fastq-gz \
	  --use-bases-mask="y100,y100" \
	  -b ${CUTOFF_LEVEL2_DIR}/post_iSAAC/inputs \
	  -o ${CUTOFF_LEVEL2_DIR}/post_iSAAC/Aligned \
	  -t ${CUTOFF_LEVEL2_DIR}/post_iSAAC/Temp \
	  -m ${MAX_RAM} &> ${CUTOFF_LEVEL2_DIR}/post_iSAAC/isaac.log

#BamCompare/puma.log: pre_iSAAC/Aligned/Projects/default/default/sorted.bam post_iSAAC/Aligned/Projects/default/default/sorted.bam
#	/illumina/development/puma/latest/bamCompare -i pre_iSAAC/Aligned/Projects/default/default/sorted.bam -i post_iSAAC/Aligned/Projects/default/default/sorted.bam -o BamCompare


simulatedVariants.vcf: EAGLE_${TEST_ID}/.L001_1.bcl.completed
	grep DEST $(dir $^)/sample_genome/canonical.vcf > $@

# "%" is for "pre" or "${CUTOFF_LEVEL2_DIR}/post"
.PRECIOUS: %_snps.vcf %_snps.compWithSim
%_snps.vcf: %_iSAAC/Aligned/Projects/default/default/sorted.bam
	${STARLING2VCF_PATH}/starlingToVcf.pl \
	  --snpsFile=$(dir $<)/Casava/Parsed/${CHROM}/snps.txt \
	  --sampleId=sample \
	  --vcfFile=$@

%_snps.compWithSim: simulatedVariants.vcf %_snps.vcf
	${EAGLE_PATH}/bin/compareVcf -s $(word 1, $^) -c $(word 2, $^) > $@ \

%_snps_aboveQ20.vcf: %_snps.vcf
	awk '{ if ($$6 >= 20) { print $0 } }' $^ > $@

%_snps_aboveQ20.compWithSim: simulatedVariants.vcf %_snps_aboveQ20.vcf
	${EAGLE_PATH}/bin/compareVcf -s $(word 1, $^) -c $(word 2, $^) > $@ \

%-comp: %_snps.compWithSim
	@echo "$@: " \
	; echo -n "Correct: " ; grep -c Correct $^ \
	; echo -n "FP     : " ; grep -c FP $^ \
	; echo -n "FN     : " ; grep -c FN $^

# "%" is for "reads" or "decodedBWT"
%.fastqWithoutReadnames.gz: %.fastq
	( sed -n 'n;p;n;n;p' $^ | gzip ) > $@

${CUTOFF_LEVEL2_DIR}/report.txt: pre_snps.compWithSim ${CUTOFF_LEVEL2_DIR}/post_snps.compWithSim reads.fastqWithoutReadnames.gz ${CUTOFF_LEVEL1_DIR}/decodedBWT.fastqWithoutReadnames.gz simulatedVariants.vcf
	@(let "CORRECT_PRE = $(shell grep -c Correct $(word 1, $^) )" \
	; let "FP_PRE = $(shell grep -c FP $(word 1, $^) )" \
	; let "FN_PRE = $(shell grep -c FN $(word 1, $^) )" \
	; let "CORRECT_POST = $(shell grep -c Correct $(word 2, $^) )" \
	; let "FP_POST = $(shell grep -c FP $(word 2, $^) )" \
	; let "FN_POST = $(shell grep -c FN $(word 2, $^) )" \
	; let "CORRECT_DIFF = CORRECT_POST - CORRECT_PRE" \
	; let "FP_DIFF = FP_POST - FP_PRE" \
	; let "FN_DIFF = FN_POST - FN_PRE" \
	; let "CORRECT_IDEAL = $(shell cut -f 2 simulatedVariants.vcf |uniq | wc -l ) - 1" \
	; let "PERCENTAGE_CORRECT_PRE = 100 * CORRECT_PRE / CORRECT_IDEAL" \
	; let "PERCENTAGE_FP_PRE = 100 * FP_PRE / CORRECT_IDEAL" \
	; let "PERCENTAGE_FN_PRE = 100 * FN_PRE / CORRECT_IDEAL" \
	; let "PERCENTAGE_CORRECT_POST = 100 * CORRECT_POST / CORRECT_IDEAL" \
	; let "PERCENTAGE_FP_POST = 100 * FP_POST / CORRECT_IDEAL" \
	; let "PERCENTAGE_FN_POST = 100 * FN_POST / CORRECT_IDEAL" \
	; let "PERCENTAGE_CORRECT_DIFF = PERCENTAGE_CORRECT_POST - PERCENTAGE_CORRECT_PRE" \
	; let "PERCENTAGE_FP_DIFF = PERCENTAGE_FP_POST - PERCENTAGE_FP_PRE" \
	; let "PERCENTAGE_FN_DIFF = PERCENTAGE_FN_POST - PERCENTAGE_FN_PRE" \
	; let "SCORE1 = PERCENTAGE_CORRECT_POST - PERCENTAGE_FP_POST - PERCENTAGE_FN_POST" \
	; let "FILESIZE_PRE = $(shell stat -L -c %s $(word 3, $^) )" \
	; let "FILESIZE_POST = $(shell stat -L -c %s $(word 4, $^) )" \
	; let "FILESIZE_RATIO = FILESIZE_POST * 100 / FILESIZE_PRE" \
	; ( echo "Final report: " \
	  ; echo " 1. Variant calls comparison \"before cutoff\" with \"after cutoff\":" \
	  ; echo "  Correct calls: $${CORRECT_PRE} -> $${CORRECT_POST} : $${CORRECT_DIFF}" \
	  ; echo "  FP           : $${FP_PRE} -> $${FP_POST} : $${FP_DIFF}" \
	  ; echo "  FN           : $${FN_PRE} -> $${FN_POST} : $${FN_DIFF}" \
	  ; echo " 2. Variant calls comparison \"before cutoff\"/\"ideal\" with \"after cutoff\"/\"ideal\":" \
	  ; echo "  Ideal correct   : $${CORRECT_IDEAL}" \
	  ; echo "  %Correct calls  : $${PERCENTAGE_CORRECT_PRE}% -> $${PERCENTAGE_CORRECT_POST}% : $${PERCENTAGE_CORRECT_DIFF}%" \
	  ; echo "  %FP             : $${PERCENTAGE_FP_PRE}% -> $${PERCENTAGE_FP_POST}% : $${PERCENTAGE_FP_DIFF}%" \
	  ; echo "  %FN             : $${PERCENTAGE_FN_PRE}% -> $${PERCENTAGE_FN_POST}% : $${PERCENTAGE_FN_DIFF}%" \
	  ; echo "  Score1 (C-FP-FN): $${SCORE1}" \
	  ; echo " 3. Comparison of bases+qualities compressibility (using gzip) \"before cutoff\" with \"after cutoff\":" \
	  ; echo "  Compressed file size: $${FILESIZE_PRE} -> $${FILESIZE_POST} : $${FILESIZE_RATIO}%" \
	  ) > $@ \
	; ( echo -e "${CUTOFF_THRESHOLD}\t${CUTOFF_QUALITY}\t$${CORRECT_PRE}\t$${FP_PRE}\t$${FN_PRE}\t$${CORRECT_POST}\t$${FP_POST}\t$${FN_POST}\t$${CORRECT_DIFF}\t$${FP_DIFF}\t$${FN_DIFF}\t$${CORRECT_IDEAL}\t$${PERCENTAGE_CORRECT_PRE}\t$${PERCENTAGE_FP_PRE}\t$${PERCENTAGE_FN_PRE}\t$${PERCENTAGE_CORRECT_POST}\t$${PERCENTAGE_FP_POST}\t$${PERCENTAGE_FN_POST}\t$${PERCENTAGE_CORRECT_DIFF}\t$${PERCENTAGE_FP_DIFF}\t$${PERCENTAGE_FN_DIFF}\t$${FILESIZE_PRE}\t$${FILESIZE_POST}\t$${FILESIZE_RATIO}" \
	  ) > $@.tsv \
	)

${CUTOFF_LEVEL2_DIR}/report_aboveQ20.txt: pre_snps_aboveQ20.compWithSim ${CUTOFF_LEVEL2_DIR}/post_snps_aboveQ20.compWithSim reads.fastqWithoutReadnames.gz ${CUTOFF_LEVEL1_DIR}/decodedBWT.fastqWithoutReadnames.gz simulatedVariants.vcf
	@(let "CORRECT_PRE = $(shell grep -c Correct $(word 1, $^) )" \
	; let "FP_PRE = $(shell grep -c FP $(word 1, $^) )" \
	; let "FN_PRE = $(shell grep -c FN $(word 1, $^) )" \
	; let "CORRECT_POST = $(shell grep -c Correct $(word 2, $^) )" \
	; let "FP_POST = $(shell grep -c FP $(word 2, $^) )" \
	; let "FN_POST = $(shell grep -c FN $(word 2, $^) )" \
	; let "CORRECT_DIFF = CORRECT_POST - CORRECT_PRE" \
	; let "FP_DIFF = FP_POST - FP_PRE" \
	; let "FN_DIFF = FN_POST - FN_PRE" \
	; let "CORRECT_IDEAL = $(shell cut -f 2 simulatedVariants.vcf |uniq | wc -l ) - 1" \
	; let "PERCENTAGE_CORRECT_PRE = 100 * CORRECT_PRE / CORRECT_IDEAL" \
	; let "PERCENTAGE_FP_PRE = 100 * FP_PRE / CORRECT_IDEAL" \
	; let "PERCENTAGE_FN_PRE = 100 * FN_PRE / CORRECT_IDEAL" \
	; let "PERCENTAGE_CORRECT_POST = 100 * CORRECT_POST / CORRECT_IDEAL" \
	; let "PERCENTAGE_FP_POST = 100 * FP_POST / CORRECT_IDEAL" \
	; let "PERCENTAGE_FN_POST = 100 * FN_POST / CORRECT_IDEAL" \
	; let "PERCENTAGE_CORRECT_DIFF = PERCENTAGE_CORRECT_POST - PERCENTAGE_CORRECT_PRE" \
	; let "PERCENTAGE_FP_DIFF = PERCENTAGE_FP_POST - PERCENTAGE_FP_PRE" \
	; let "PERCENTAGE_FN_DIFF = PERCENTAGE_FN_POST - PERCENTAGE_FN_PRE" \
	; let "SCORE1 = PERCENTAGE_CORRECT_POST - PERCENTAGE_FP_POST - PERCENTAGE_FN_POST" \
	; let "FILESIZE_PRE = $(shell stat -L -c %s $(word 3, $^) )" \
	; let "FILESIZE_POST = $(shell stat -L -c %s $(word 4, $^) )" \
	; let "FILESIZE_RATIO = FILESIZE_POST * 100 / FILESIZE_PRE" \
	; ( echo "Final report: " \
	  ; echo " 1. Variant calls comparison \"before cutoff\" with \"after cutoff\":" \
	  ; echo "  Correct calls: $${CORRECT_PRE} -> $${CORRECT_POST} : $${CORRECT_DIFF}" \
	  ; echo "  FP           : $${FP_PRE} -> $${FP_POST} : $${FP_DIFF}" \
	  ; echo "  FN           : $${FN_PRE} -> $${FN_POST} : $${FN_DIFF}" \
	  ; echo " 2. Variant calls comparison \"before cutoff\"/\"ideal\" with \"after cutoff\"/\"ideal\":" \
	  ; echo "  Ideal correct   : $${CORRECT_IDEAL}" \
	  ; echo "  %Correct calls  : $${PERCENTAGE_CORRECT_PRE}% -> $${PERCENTAGE_CORRECT_POST}% : $${PERCENTAGE_CORRECT_DIFF}%" \
	  ; echo "  %FP             : $${PERCENTAGE_FP_PRE}% -> $${PERCENTAGE_FP_POST}% : $${PERCENTAGE_FP_DIFF}%" \
	  ; echo "  %FN             : $${PERCENTAGE_FN_PRE}% -> $${PERCENTAGE_FN_POST}% : $${PERCENTAGE_FN_DIFF}%" \
	  ; echo "  Score1 (C-FP-FN): $${SCORE1}" \
	  ; echo " 3. Comparison of bases+qualities compressibility (using gzip) \"before cutoff\" with \"after cutoff\":" \
	  ; echo "  Compressed file size: $${FILESIZE_PRE} -> $${FILESIZE_POST} : $${FILESIZE_RATIO}%" \
	  ) > $@ \
	; ( echo -e "${CUTOFF_THRESHOLD}\t${CUTOFF_QUALITY}\t$${CORRECT_PRE}\t$${FP_PRE}\t$${FN_PRE}\t$${CORRECT_POST}\t$${FP_POST}\t$${FN_POST}\t$${CORRECT_DIFF}\t$${FP_DIFF}\t$${FN_DIFF}\t$${CORRECT_IDEAL}\t$${PERCENTAGE_CORRECT_PRE}\t$${PERCENTAGE_FP_PRE}\t$${PERCENTAGE_FN_PRE}\t$${PERCENTAGE_CORRECT_POST}\t$${PERCENTAGE_FP_POST}\t$${PERCENTAGE_FN_POST}\t$${PERCENTAGE_CORRECT_DIFF}\t$${PERCENTAGE_FP_DIFF}\t$${PERCENTAGE_FN_DIFF}\t$${FILESIZE_PRE}\t$${FILESIZE_POST}\t$${FILESIZE_RATIO}" \
	  ) > $@.tsv \
	)



clean:
	rm -rf "EAGLE_${TEST_ID}" FastQ reads.fa reads.fastq BWT outFileEndPos.bwt decodedBWT.fa* pre_iSAAC post_iSAAC pre_snps.txt post_snps.txt
