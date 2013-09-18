BEETL: Burrows-Wheeler Extended Tool Library
============================================

Copyright (c) 2013 Illumina, Inc.

This software is covered by the "Illumina Non-Commercial Use Software
and Source Code License Agreement" and any user of this software or
source file is bound by the terms therein.


Description
-----------

BEETL is a suite of applications for building and manipulating the 
Burrows-Wheeler Transform (BWT) of collections of DNA sequences.  
The algorithms employed in BEETL are intended to scale to collections of 
sequences containing one billion entries or more.


Installation
------------

Assuming the following paths:
- /sourcePath : BEETL source code directory
- /buildPath : temporary build directory 
- /installPath : final installation directory

If you obtained the source from git:

    cd /sourcePath
    ./bootstrap_with_autotools

If you have gcc>=4.7, you can activate OpenMP parallelism with:

    export CXXFLAGS="-DUSE_OPENMP -fopenmp"

Then, for everyone:

    cd /buildPath
    /sourcePath/configure --prefix=/installPath
    make
    make install
    export PATH=$PATH:/installPath/bin


Tools
-----

Running `beetl --help` shows a high-level description of each BEETL tool:

    beetl bwt       Encode a set of nucleotide sequences using Burrows-Wheeler transform
    beetl unbwt     Decode nucleotide sequences dataset from its BWT
    beetl correct   Correct sequencing errors in a BWT
    beetl compare   Compare two BWT datasets
    beetl search    Search within a BWT dataset
    beetl convert   Convert between file formats


Examples
--------

### FASTA to BWT to FASTA

    # FASTA to BWT
    beetl-bwt -i input.fasta -o myBWT --output-format=ascii
    
    # BWT to FASTA
    beetl-unbwt -i myBWT -o output.fasta


### Error correction using BWT

    # BWT creation from FASTA
    beetl-bwt -i input.fasta --add-rev-comp --generate-end-pos-file
    
    # BWT correction, generating a list of corrections
    beetl-correct -i outBwt -o corrections.csv -L 2000000 -e 10000 -k 30 -w 13

    # Applying corrections to origing fasta file
    align-corrector-strings -i input.fasta -c corrections.csv -o corrected.fasta --input-reads-format=fasta --output-reads-format=fasta -a no-indels -q '?' --min-witness-length=14


### Tumour-normal comparison using BWT

    # Generation of tumour and normal BWTs
    beetl-bwt -i tumour.fastq -o tumourBWT --generate-end-pos-file
    beetl-bwt -i normal.fastq -o normalBWT --generate-end-pos-file
    
    # BWT comparison with default output to stdout
    beetl-compare --mode=split -a tumourBWT -b normalBWT


### Meta-BEETL Metagenomics

We describe here how to reproduce the results presented in our paper [metaBEETL: high-throughput analysis of heterogeneous microbial populations from shotgun DNA sequences](http://www.biomedcentral.com/1471-2105/14/S5/S2).  
The first challenge is to get the metagenomic database. You can either download it or rebuild it yourself:

#### Downloading the metagenomic database

Note: If you are at Illumina, you can use the pre-installed data: `export METAGENOME_DATABASE_PATH=/illumina/scratch/BWT/metagenomics/allNucleotides/workDirFullRleInt`

Otherwise, all the references and associated metadata used in the paper are available from Amazon S3:

    wget https://s3.amazonaws.com/metaBEETL/metaBEETL.ncbi.refseq.tar.gz    # 66 GB

After downloading this file, decompress it in a directory which we will then refer to as `${METAGENOME_DATABASE_PATH}`

    mkdir BeetlMetagenomeDatabase
    cd BeetlMetagenomeDatabase
    tar xzf ../metaBEETL.ncbi.refseq.tar.gz
    export METAGENOME_DATABASE_PATH=`pwd`

If you downloaded those files, you can skip the next section.


#### Building the metagenomic database

How to create a database of reference genomes for metaBEETL: 

1. Creating a new metagenomic database requires an installation of the SeqAn library (www.seqan.de).  
   You can set the location of your SeqAn installation before compiling BEETL by using the configure parameter --with-seqan.  
   If you do this, all executables mentioned below will be compiled automatically and copied into /installPath/bin.

2. Download all bacteria, archea and virus genomes from ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/.

3. Create single sequence files with the reverse complement of all genomes but exclude plasmids.  
   For this, use genomeToSingleSeq.  
   Do NOT change the names of the generated files.

4. Open the script arrayBWT.sh and adapt the paths.  
   To create the BWTs for all the G\_* and G\_*\_rev files using a Grid Engine cluster, submit `qsub -t n arrayBWT.sh`, where n should be 1-500 if you have the files G\_1 till G\_500.  
   As an alternative you can also run makeBWTSkew for each of the files individually.

5. Login on a machine with enough memory to load all sequences in RAM and run mergeBacteria on all files.  
   You will have to call this once for each of the piles created by BEETL which means six times overall e.g.  
```
   for pileNum in `seq 0 5`; do mergeBacteria $pileNum databaseTest G_* ; done
```  
   For each pile this will create 3 files and one fileCounter.csv (which doesn't change)  
   The output files will be called: yourPrefix-A0\*, -B0\* and -C0*  
   - C0* contain, for each position in the BWT, the file number where the char at this position came from.  
   - B0* are the BWTs for all files.  
   - A0* contain, for each position in the BWT, the suffix position in the file where this BWT came from.  
  
   Files C0* and B0* are needed for the countWords algorithm. 

6. Download the NCBI taxonomy from here. ftp://ftp.ncbi.nih.gov/pub/taxonomy/.  
   You will need the files names.dmp, nodes.dmp, merged.dmp and the gi\_taxid\_nucl.cmp.  
   Use the findTaxa script to find the taxonomic tree corresponding to the file numbers in the database.  
   You will need the headerfile produced by running "genomeToSingleSeq" and 
    fileCounterFile created during the merging of the bacterial reference genomes.  
    This will first give all found GI-Numbers out in std::err.  
    If you want to shorten the search for the gi numbers in gi\_taxid\_nucl.dmp, grep the gi-numbers in advance and search for 
    the taxonomy only in the grepped output.  
    Finally, you get for each file number in the database a taxonomic tree with the taxonomic ids.  
    There will be some 0 in the taxonomic tree. This is a taxonomic id which could not be 
    matched to: Superkingdom, Phylum, Order, Class, Family, Genus, Species or Strain.  
    Sometimes there are just missing taxa in the taxonomy. 

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

Create the BWT of the input dataset, with reverse-complemented reads:

    beetl-bwt -i SRS013948.seq -o bwt_SRS013948 --add-rev-comp

Run BEETL metagenomic classification:

    beetl-compare \
      --mode=metagenomics \
      -a bwt_SRS013948 \
      -b ${METAGENOME_DATABASE_PATH}/ncbiMicros \
      -t ${METAGENOME_DATABASE_PATH}/ncbiFileNumToTaxTree \
      -w 50 \
      -n 1 \
      -k 100 \
        > metaBeetlOutput.txt

Since the algorithm repeatedly looks up the filenumbers for each BWT-Position we recommend to put these files on a disk with fast read access.  
Setting k = sequence length gets you the maximal amount of information but the output file will be large.


#### Visualisation of metagenomic results

The `parseMetagenomeOutput` tool is used to obtain a summary of the metaBEETL output which can then be visualized.

Run:

    parseMetagenomeOutput \
      -f \
      -b metaBeetlOutput.txt \
      -t ${METAGENOME_DATABASE_PATH}/ncbiFileNumToTaxTree \
      -w 50 75 100

-w is a vector of wordSizes for which results should be generated.  
The results of the parsing are saved in files named after the chosen word sizes.  
In addition a file will be generated with the information about the largest BWT positions.

You can look at the files `50`, `75` and `100` by hand, and if you are at Illumina, you can
visualise them by running the java program:

     /illumina/scratch/BWT/metagenomics/Visualisation/visualisation

Click on the bottom-left button "Load results" to open a data file.

