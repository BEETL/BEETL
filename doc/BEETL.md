BEETL: Burrows-Wheeler Extended Tool Library
============================================

Copyright (c) 2011-2014 Illumina, Inc.

This software is covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
and any user of this software or source file is bound by the terms therein.


Description
-----------

BEETL is a suite of applications for building and manipulating the 
Burrows-Wheeler Transform (BWT) of collections of DNA sequences.  
The algorithms employed in BEETL are intended to scale to collections of 
sequences containing one billion entries or more.


External dependencies
---------------------

On a clean Amazon Ubuntu instance, you would need:

    sudo apt-get update
    sudo apt-get install git unzip g++ make automake


On Mac with LLVM-clang (although unsupported as it doesn't have OpenMP), you may need to define `export CXX=c++; export CXXFLAGS="-stdlib=libc++"`

On Mac with g++ (try macports) installed in /opt/local (warning: it may be in /usr/local), you may need to define (adjust for your gcc version and path) `export LDFLAGS="-L/opt/local/lib/gcc49"; export CPPFLAGS="-I/opt/local/include/gcc49"`


Installation
------------

Assuming the following paths:
- /sourcePath : BEETL source code directory
- /buildPath : temporary build directory 
- /installPath : final installation directory

```sh
    cd /buildPath
    /sourcePath/configure --prefix=/installPath
    make
    make install
    export PATH=$PATH:/installPath/bin
```

Tools
-----

Running `beetl --help` shows a high-level description of each BEETL tool:

    beetl bwt       Encode a set of nucleotide sequences using Burrows-Wheeler transform
    beetl unbwt     Decode nucleotide sequences dataset from its BWT
    beetl correct   Correct sequencing errors in a BWT
    beetl compare   Compare two BWT datasets
    beetl search    Search within a BWT dataset
    beetl extend    Extend BWT intervals to identify their associated sequence numbers
    beetl index     Generate index for BWT file to speed other algorithms up
    beetl convert   Convert between file formats


...and some higher level tools:

    beetl-fastq     Index and search FASTQ files
    beetl-flow-tumour-normal-fastq-filter   Tumour-normal subtraction, filtering out uninteresting reads for faster downstream analysis


For computing LCP array
-------------------------------------
Set BUILD_LCP to 1 in src/shared/Tools.hh and compile

Example:

    beetl-bwt -i input.fasta -o myBWT --output-format=ascii --concatenate-output --generate-lcp
    

Paired and reverse-complemented reads
-------------------------------------

If using paired reads, two input styles can be used:
- the two reads can be provided in an "all 1 all 2" configuration (rather than interleaved): all reads 1 as a first set, and all reads 2 as a second set, concatenated in one input file. `beetl-bwt --paired-reads-input` can then be used to keep track of this configuration.
- or each two paired reads can be concatenated as one and `beetl-bwt --sub-sequence-length` will split them.

Reverse-complemented reads are not usually present in the input file, but generated during BWT construction by `beetl-bwt --add-rev-comp`.


Examples
--------

### FASTA to BWT to FASTA

    # FASTA to BWT
    beetl-bwt -i input.fasta -o myBWT --output-format=ascii
    
    # BWT to FASTA
    beetl-unbwt -i myBWT -o output.fasta


### K-mer search in FastQ using BWT (bases only, see below for qualities)

    # FastQ to BWT
    beetl-bwt -i input.fastq -o bwt --generate-end-pos-file

    # Optional: Index BWT files
    beetl-index -i bwt

    # K-mer search
    beetl-search -i bwt -k ACGT -o searchedKmers.bwtIntervals

    # Propagation of the found intervals to end of reads, in order to identify sequence numbers
    beetl-extend -i searchedKmers.bwtIntervals -b bwt -o searchedKmers.sequenceNumbers

    # Extraction of the FastQ lines
    beetl-convert -i input.fastq --extract-sequences=searchedKmers.sequenceNumbers -o sequencesWithSearchedKmers.fastq


### K-mer search in FastQ (returning full reads: read ids + bases + quality scores)

    # FastQ conversion
    beetl-fastq --mode=init -i myFile.fastq -o myIndexedFile

    # K-mer search + extraction of the FastQ lines
    beetl-fastq --mode=search -i myIndexedFile -k ACGT -o ACGT.fastq


### Error correction using BWT

    # BWT creation from FASTA
    beetl-bwt -i input.fasta --add-rev-comp --generate-end-pos-file
    
    # BWT correction, generating a list of corrections
    beetl-correct -i outBwt -o corrections.csv -L 2000000 -e 10000 -k 30 -w 13

    # Applying corrections to origing fasta file
    beetl-correct-apply-corrections -i input.fasta -c corrections.csv -o corrected.fasta --input-reads-format=fasta --output-reads-format=fasta -a no-indels -q '?' --min-witness-length=14


### Tumour-normal filtering using BWT

The tool 'beetl-flow-tumour-normal-fastq-filter' takes 2 datasets (tumour and normal), each made of 2 paired-end fastq or fastq.gz files (one for read 1, one for read 2).  
It extracts all the read pairs that contain suspected variants, and outputs 4 fastq files, subsets of the input files.  
The output files are placed in the specified output directory and are called:
- tumour\_read1.fastq
- tumour\_read2.fastq
- normal\_read1.fastq
- normal\_read2.fastq

Example:

    beetl-flow-tumour-normal-fastq-filter \
      Tumour/lane1_NoIndex_L001_R1_001.fastq.gz \
      Tumour/lane1_NoIndex_L001_R2_001.fastq.gz \
      Normal/lane1_NoIndex_L001_R1_001.fastq.gz \
      Normal/lane1_NoIndex_L001_R2_001.fastq.gz \
      filteredOutputDirectory

Note: Intermediate files (also present in the specified output directory) are not currently deleted, as they may be useful for further experiments.


### Meta-BEETL Metagenomics

We describe here how to reproduce the results presented in our paper [metaBEETL: high-throughput analysis of heterogeneous microbial populations from shotgun DNA sequences](http://www.biomedcentral.com/1471-2105/14/S5/S2).  
The first challenge is to get the metagenomic database. You can either download it or rebuild it yourself:

#### Downloading the metagenomic database

Note: If you are at Illumina, you can use the pre-installed data: `export METAGENOME_DATABASE_PATH=/illumina/scratch/BWT/metagenomics/metaBeetlNcbiDb/stable` (Note: use the `testing` version if you are also using the testing version of the tools)

Otherwise, all the references and associated metadata used in the paper are available from Amazon S3 (24GB of files):

    mkdir BeetlMetagenomeDatabase
    cd BeetlMetagenomeDatabase
    for i in `curl https://s3.amazonaws.com/metaBEETL | grep -oP "[^>]*bz2"` ; \
        do wget https://s3.amazonaws.com/metaBEETL/$i & done

After downloading these files, decompress them in a directory which we will then refer to as `${METAGENOME_DATABASE_PATH}`

    bunzip2 *.bz2
    export METAGENOME_DATABASE_PATH=`pwd`

If you downloaded those files, you can skip the next section.


#### Building the metagenomic database

How to create a database of reference genomes for metaBEETL: 

1. Creating a new metagenomic database requires an installation of the SeqAn library (www.seqan.de).  
   You can set the location of your SeqAn installation before compiling BEETL by using the configure parameter --with-seqan.  
   If you do this, all executables mentioned below will be compiled automatically and copied into /installPath/bin.

2. Download all bacteria, archea and virus genomes from ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/
```
   mkdir metaBeetlDb
   cd metaBeetlDb
   wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/all.fna.tar.gz -P downloads      # 2.7 GB
   
   mkdir all.fna
   tar xzf downloads/all.fna.tar.gz -C all.fna/
```

3. Create single sequence files with the reverse complement of all genomes but exclude plasmids.  
   For this, use metabeetl-db-genomesToSingleSeq.  
   Do NOT change the names of the generated files.
```
   metabeetl-db-genomesToSingleSeq -f all.fna/*/*.fna -s singleSeqGenomes
   cd singleSeqGenomes
```

4. Open the script metabeetl-db-arrayBWT.sh and adapt the paths.  
   To create the BWTs for all the G\_\* and G\_\*\_rev files using a Grid Engine cluster, submit `qsub -t n metabeetl-db-arrayBWT.sh`, where n should be 1-500 if you have the files G\_1 till G\_500.  
   As an alternative you can also run metabeetl-db-makeBWTSkew for each of the files individually.
```
   cp `which metabeetl-db-arrayBWT.sh` .
   vi metabeetl-db-arrayBWT.sh  # Adjust paths in this file
   qsub -t n metabeetl-db-arrayBWT.sh
```
 OR
```
   ( echo -n "all: " ; for i in G_*; do echo -n " bwt_${i}-B00"; done ; echo -e "\n" ; \
      for i in G_*; do echo "bwt_${i}-B00: ${i}"; echo -e "\tmetabeetl-db-makeBWTSkew ${i} ${i}\n" ; done ) > Makefile
   make -j
```

5. Login on a machine with enough memory to load all sequences in RAM (~60GB) and run mergeBacteria on all files.  
   You will have to call this once for each of the piles created by BEETL which means six times overall e.g.  
```
   for pileNum in `seq 0 5`; do metabeetl-db-mergeBacteria $pileNum ncbiMicros <( ls G_* ) ; done
```  
   For each pile this will create 3 files and one fileCounter.csv (which doesn't change)  
   The output files will be called: ncbiMicros-A0\*, -B0\* and -C0\*  
   - -A0\* contain, for each position in the BWT, the suffix position in the file where this BWT came from.  
   - -B0\* are the BWTs for all files.  
   - -C0\* contain, for each position in the BWT, the file number where the char at this position came from.  
  
   Files -B0\* and -C0\* are needed for the countWords algorithm. 

6. Optionally, convert the ASCII BWT files to RLE BWT, to make subsequent processing faster:  
```
   for pileNum in `seq 0 5`; do \  
     mv ncbiMicros-B0${pileNum} ncbiMicros-B0${pileNum}.ascii ;  \
     beetl-convert \  
       --input-format=bwt_ascii \  
       --output-format=bwt_rle \  
       -i ncbiMicros-B0${pileNum}.ascii \  
       -o ncbiMicros-B0${pileNum} ; \  
   done
```

7. Download the NCBI taxonomy from ftp://ftp.ncbi.nih.gov/pub/taxonomy/  
   You will need the files names.dmp, nodes.dmp and the gi\_taxid\_nucl.dmp, which are packaged as taxdump.tar.gz and gi\_taxid\_nucl.dmp.gz:  
```
   cd downloads
   wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz         # 26 MB
   wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz   # 800 MB

   tar xzf taxdump.tar.gz names.dmp nodes.dmp
   gunzip gi_taxid_nucl.dmp.gz
   cd ..
```
   Use the metabeetl-db-findTaxa script to find the taxonomic tree corresponding to the file numbers in the database.  
   You will need the headerFile produced by running "metabeetl-db-genomesToSingleSeq" and 
   fileCounter created during the merging of the bacterial reference genomes.  
   Finally, you get for each file number in the database a taxonomic tree with the taxonomic ids.  
   There will be some 0 in the taxonomic tree. This is a taxonomic id which could not be 
   matched to: Superkingdom, Phylum, Order, Class, Family, Genus, Species or Strain.  
   Sometimes there are just missing taxa in the taxonomy. We supplement this with the file `metaBeetlExtraNames.dmp` below. 
```
   metabeetl-db-findTaxa \
     -nA downloads/names.dmp \
     -nO downloads/nodes.dmp \
     -nG downloads/gi_taxid_nucl.dmp \
     -h singleSeqGenomes/headerFile.csv \
     -f singleSeqGenomes/filecounter.csv \
     > ncbiFileNumToTaxTree

   ( grep scientific downloads/names.dmp ; cat ${BEETL_INSTALL_DIR}/share/beetl/metaBeetlExtraNames.dmp ) > metaBeetlTaxonomyNames.dmp
```

8. (Experimental, but nice) Run meta-BEETL on each of the genomes to calculate the normalisation factors
```
   mkdir normalisation
   cd normalisation
   for genome in ../singleSeqGenomes/G_*; do
     (
      genomeNum=`basename ${genome}`
      mkdir ${genomeNum}
      cd ${genomeNum}
      for i in ../../singleSeqGenomes/bwt_${genomeNum}-B0?; do ln -s ${i} ; done
      for i in `seq 0 6`; do touch bwt_${genomeNum}-B0${i}; done
      time beetl-compare --mode=metagenomics -a bwt_${genomeNum} -b ../../singleSeqGenomes/ncbiMicros -t ../../ncbiFileNumToTaxTree -w 20 -n 1 -k 50 --no-comparison-skip -v &> out.step1
      rm -f BeetlCompareOutput/cycle51.subset*
                  touch empty.txt
      time cat BeetlCompareOutput/cycle*.subset* | metabeetl-convertMetagenomicRangesToTaxa ../../ncbiFileNumToTaxTree ../../singleSeqGenomes/ncbiMicros ../../metaBeetlTaxonomyNames.dmp empty.txt 20 50 - &> out.step2
      cd ..
     ) &
   done ; wait

   for i in `seq 1 2977`; do echo "G_$i"; X=`grep -P "G_${i}$" ../singleSeqGenomes/filecounter.csv |cut -f 1 -d ','`; TAX=`grep -P "^${X} " ../ncbiFileNumToTaxTree | tr -d "\n\r"` ; echo "TAX(${X}): ${TAX}."; TAXID=`echo "${TAX}" | sed 's/\( 0\)*$//g' |awk '{print $NF}'`; echo "TAXID=${TAXID}"; COUNTS=`grep Counts G_${i}/out.step2 | head -1`; echo "COUNTS=${COUNTS}"; MAIN_COUNT=`echo "${COUNTS}  " | sed "s/^.* ${TAXID}:\([0-9]*\) .*$/\1/ ; s/Counts.*/0/"` ; echo "MAIN_COUNT=${MAIN_COUNT}" ; SUM=`echo "${COUNTS}  " | tr ' ' '\n' | sed 's/.*://' | awk 'BEGIN { sum=0 } { sum+=$1 } END { print sum }'` ; echo "SUM=$SUM"; PERCENT=`echo -e "scale=5\n100*${MAIN_COUNT}/${SUM}" | bc` ; echo "PERCENT=${PERCENT}" ; echo "FINAL: G_${i} ${TAXID} ${MAIN_COUNT} ${SUM} ${COUNTS}" ; done > r2977

   grep FINAL r2977 > ../normalisation.txt
```

9. Link (or move) all the useful files to a final directory
```
   mkdir metaBeetlNcbiDb
   cd metaBeetlNcbiDb
   ln -s ../ncbiFileNumToTaxTree
   ln -s ../normalisation.txt
   ln -s ../downloads/metaBeetlTaxonomyNames.dmp
   ln -s ../singleSeqGenomes/filecounter.csv
   ln -s ../singleSeqGenomes/headerFile.csv
   for i in ../singleSeqGenomes/ncbiMicros-[BC]0[0-5]; do ln -s $i ; done
```

Now everything should be ready to run metaBEETL!  :-)

#### Downloading the metagenomic input dataset

The sample SRS013948 from the Human Microbiome project should be available from:

http://downloads.hmpdacc.org/data/Illumina/throat/SRS013948.tar.bz2  
or  
ftp://public-ftp.hmpdacc.org/Illumina/throat/SRS013948.tar.bz2

Decompressing this file

    tar xjf SRS013948.tar.bz2
    
should give you:

    SRS013948.denovo_duplicates_marked.trimmed.1.fastq
    SRS013948.denovo_duplicates_marked.trimmed.2.fastq
    SRS013948.denovo_duplicates_marked.trimmed.singleton.fastq

We need to normalise the lengths of these sequences and put them all in one file:

    beetl-convert \
      -i SRS013948.denovo_duplicates_marked.trimmed.1.fastq \
      -o paddedSeq1.seq \
      --sequence-length=100
    beetl-convert \
      -i SRS013948.denovo_duplicates_marked.trimmed.2.fastq \
      -o paddedSeq2.seq \
      --sequence-length=100
    beetl-convert \
      -i SRS013948.denovo_duplicates_marked.trimmed.singleton.fastq \
      -o paddedSeqSingleton.seq \
      --sequence-length=100
    cat paddedSeq1.seq paddedSeq2.seq paddedSeqSingleton.seq > SRS013948.seq


#### Running BEETL in metagenomic mode

Create the BWT of the input dataset:

    beetl-bwt -i SRS013948.seq -o bwt_SRS013948

Run BEETL metagenomic classification:

    beetl-compare \
      --mode=metagenomics \
      -a bwt_SRS013948 \
      -b ${METAGENOME_DATABASE_PATH}/ncbiMicros \
      -t ${METAGENOME_DATABASE_PATH}/ncbiFileNumToTaxTree \
      -w 20 \
      -n 1 \
      -k 50 \
      --no-comparison-skip

It currently generates many output files in a `BeetlCompareOutput` directory, organised per cycle and containing the information about the BWT ranges of k-mers that start diverging between the input dataset and the database.
Setting k = sequence length gets you the maximal amount of information but the output file will be large.


#### Gathering and visualisation of metagenomic results

The `metabeetl-convertMetagenomicRangesToTaxa` tool converts the BWT ranges of k-mer matches from the previous stage into genome and ancestor IDs and generates text and graphical output files (currently html files using the Krona javascript visualisation library).

Since the algorithm repeatedly looks up the filenumbers for each BWT position we recommend to put these ncbiMicros-C0* files on a disk with fast read access.
Alternatively, if you have enough RAM, you can try `metabeetl-convertMetagenomicRangesToTaxa_withMmap`, which maps these files to RAM.

Run:

    cat BeetlCompareOutput/cycle*.subset* | \
    metabeetl-convertMetagenomicRangesToTaxa \
      ${METAGENOME_DATABASE_PATH}/ncbiFileNumToTaxTree \
      ${METAGENOME_DATABASE_PATH}/ncbiMicros \
      ${METAGENOME_DATABASE_PATH}/metaBeetlTaxonomyNames.dmp \
      ${METAGENOME_DATABASE_PATH}/normalisation.txt \
      20 \
      50 \
      - \
      > metaBeetl.log

Three TSV (tab-separated values) files are generated (column 1: Taxonomy Id, column 2: Taxonomy level, column 3: k-mer count, column 4: k-mer count including children):
- metaBeetl.tsv: raw k-mer counts for every leaves and ancestors.
- metaBeetl_normalised.tsv: some counts from ancestors are moved towards leaf items. The proportion thereof is pre-computed by aligning each individual genome to the full database.
- metaBeetl_normalised2.tsv: Only leaves
 of the taxonomy tree are kept, and counts are normalised relatively to genome sizes. 

Three Krona files are also generated:
- metaBeetl_krona.html
- metaBeetl_krona_normalised.html
- metaBeetl_krona_normalised2.html
