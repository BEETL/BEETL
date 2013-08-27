SUBDIRS = . frontends metagenomics

bin_PROGRAMS = OldBeetl
noinst_LIBRARIES = libzoo.a liball.a

OldBeetl_SOURCES = shared/Beetl.cpp
OldBeetl_LDADD = liball.a libzoo.a ${BOOST_LDADD}
OldBeetl_LDFLAGS = ${BOOST_LDFLAGS}
OldBeetl_CXXFLAGS = -I$(srcdir) -I$(srcdir)/shared -I$(srcdir)/BCR -I$(srcdir)/BCRext -I$(srcdir)/countWords

# Unfinished tool?
#indexBwt_SOURCES = BCRext/IndexBwt.cpp
#indexBwt_LDADD = liball.a
#indexBwt_CXXFLAGS = -I$(srcdir) -I$(srcdir)/shared -I$(srcdir)/BCR -I$(srcdir)/BCRext -I$(srcdir)/countWords


liball_a_SOURCES = \
	BCRext/BCRext.cpp \
	BCR/BCR_BWTCollection.cpp \
	countWords/CountWords.cpp \
	shared/Timer.cpp \
	shared/Tools.cpp \
	BCR/BWTCollection.cpp \
	BCR/BCRexternalBWT.cpp \
	BCR/BuildBCR.cpp \
	BCR/PredictiveEncoding.cpp \
	BCR/Sorting.cpp \
	BCR/TransposeFasta.cpp \
	BCRext/BwtReader.cpp \
	BCRext/BwtWriter.cpp \
	BCRext/ReadBuffer.cpp \
	shared/Filename.cpp \
	shared/SeqReader.cpp \
	countWords/BackTracker.cpp \
	countWords/IntervalHandlerBase.cpp \
	countWords/IntervalHandlerMetagenome.cpp \
	countWords/IntervalHandlerReference.cpp \
	countWords/IntervalHandlerSplice.cpp \
	countWords/RangeStore.cpp

liball_a_CXXFLAGS = -I$(srcdir) -I$(srcdir)/shared -I$(srcdir)/BCR -I$(srcdir)/BCRext -I$(srcdir)/countWords ${BOOST_CXXFLAGS}


libzoo_a_SOURCES = \
	libzoo/cli/Common.cpp \
	libzoo/cli/ToolParameters.cpp \
	libzoo/io/Bcl.cpp \
	libzoo/io/FastOFStream.cpp \
	libzoo/util/Logger.cpp \
	libzoo/util/TemporaryFilesManager.cpp

libzoo_a_CXXFLAGS = -I$(srcdir) ${BOOST_CXXFLAGS}


# Unused file?
# Transpose.o: BCR/Transpose.cpp
