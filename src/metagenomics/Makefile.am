bin_PROGRAMS = metabeetl-db-genomesToSingleSeq metabeetl-db-makeBWTSkew metabeetl-db-mergeBacteria metabeetl-db-findTaxa metabeetl-parseMetagenomeOutput metabeetl-convertMetagenomicRangesToTaxa metabeetl-convertMetagenomicRangesToTaxa_withMmap
dist_bin_SCRIPTS = metabeetl-db-arrayBWT.sh


metabeetl_db_genomesToSingleSeq_SOURCES = GenomesToSingleSeq.cpp

metabeetl_db_makeBWTSkew_SOURCES = BuildChromosomeBwt.cpp
metabeetl_db_makeBWTSkew_CXXFLAGS = ${SEQAN_CXXFLAGS}

metabeetl_db_mergeBacteria_SOURCES = MergeBacteria.cpp
metabeetl_db_mergeBacteria_CXXFLAGS = -I$(srcdir) -I$(top_srcdir)/src/shared 

metabeetl_db_findTaxa_SOURCES = findCertainTaxLevel.cpp
metabeetl_db_findTaxa_LDADD = ../liball.a ../libzoo.a

metabeetl_convertMetagenomicRangesToTaxa_SOURCES = ConvertMetagenomicRangesToTaxa.cpp
metabeetl_convertMetagenomicRangesToTaxa_CXXFLAGS = ${OPENMP_CXXFLAGS}
metabeetl_convertMetagenomicRangesToTaxa_LDADD = ../liball.a ../libzoo.a

metabeetl_convertMetagenomicRangesToTaxa_withMmap_SOURCES = ConvertMetagenomicRangesToTaxa.cpp
metabeetl_convertMetagenomicRangesToTaxa_withMmap_CXXFLAGS = ${OPENMP_CXXFLAGS} -DUSE_MMAP
metabeetl_convertMetagenomicRangesToTaxa_withMmap_LDADD = ../liball.a ../libzoo.a

# Soon to be deprecated:
metabeetl_parseMetagenomeOutput_SOURCES = parse.cpp
metabeetl_parseMetagenomeOutput_LDADD = ../liball.a ../libzoo.a


noinst_HEADERS = metaShared.hh Krona.hh OutputTsv.hh
