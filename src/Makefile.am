SUBDIRS = . frontends metagenomics

bin_PROGRAMS = OldBeetl
noinst_LIBRARIES = liball.a

OldBeetl_SOURCES = shared/Beetl.cpp
OldBeetl_LDADD = liball.a
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
	BCR/Sorting.cpp \
	BCR/TransposeFasta.cpp \
	BCRext/BwtReader.cpp \
	BCRext/BwtWriter.cpp \
	BCRext/ReadBuffer.cpp \
	shared/SeqReader.cpp \
	shared/TemporaryFilesManager.cpp \
	shared/Logger.cpp \
	countWords/BackTracker.cpp \
	countWords/IntervalHandlerBase.cpp \
	countWords/IntervalHandlerMetagenome.cpp \
	countWords/IntervalHandlerReference.cpp \
	countWords/IntervalHandlerSplice.cpp \
	countWords/RangeStore.cpp

liball_a_CXXFLAGS = -I$(srcdir) -I$(srcdir)/shared -I$(srcdir)/BCR -I$(srcdir)/BCRext -I$(srcdir)/countWords


# Unused file?
# Transpose.o: BCR/Transpose.cpp
