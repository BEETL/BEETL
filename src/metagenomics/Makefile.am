bin_PROGRAMS = genomesToSingleSeq makeBWTSkew mergeBacteria findTaxa parseMetagenomeOutput convertMetagenomicRangesToTaxa
dist_bin_SCRIPTS = arrayBWT.sh


genomesToSingleSeq_SOURCES = GenomesToSingleSeq.cpp

makeBWTSkew_SOURCES = BuildChromosomeBwt.cpp
makeBWTSkew_CXXFLAGS = ${SEQAN_CXXFLAGS}

mergeBacteria_SOURCES = MergeBacteria.cpp
mergeBacteria_CXXFLAGS = -I$(srcdir) -I$(top_srcdir)/src/shared 

findTaxa_SOURCES = findCertainTaxLevel.cpp
findTaxa_LDADD = ../liball.a ../libzoo.a

parseMetagenomeOutput_SOURCES = parse.cpp
parseMetagenomeOutput_LDADD = ../liball.a ../libzoo.a

convertMetagenomicRangesToTaxa_SOURCES = ConvertMetagenomicRangesToTaxa.cpp
convertMetagenomicRangesToTaxa_LDADD = ../liball.a ../libzoo.a


noinst_HEADERS = metaShared.hh Krona.hh
