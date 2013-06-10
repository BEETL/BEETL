bin_PROGRAMS = genomesToSingleSeq makeBWTSkew mergeBacteria findTaxa parseMetagenomeOutput
bin_SCRIPTS = arrayBWT.sh


genomesToSingleSeq_SOURCES = GenomesToSingleSeq.cpp

makeBWTSkew_SOURCES = BuildChromosomeBwt.cpp
makeBWTSkew_CXXFLAGS = ${SEQAN_CXXFLAGS}

mergeBacteria_SOURCES = MergeBacteria.cpp
mergeBacteria_CXXFLAGS = -I$(srcdir) -I$(top_srcdir)/src/shared 


findTaxa_SOURCES = findCertainTaxLevel.cpp

parseMetagenomeOutput_SOURCES = parse.cpp
