SUBDIRS = . frontends metagenomics

AM_CXXFLAGS = ${OPENMP_CXXFLAGS}

bin_PROGRAMS = OldBeetl
noinst_LIBRARIES = libzoo.a liball.a

OldBeetl_SOURCES = shared/Beetl.cpp shared/Beetl.hh
OldBeetl_LDADD = liball.a libzoo.a ${BOOST_LDADD}
OldBeetl_LDFLAGS = ${BOOST_LDFLAGS}
OldBeetl_CXXFLAGS = -I$(srcdir) -I$(srcdir)/shared -I$(srcdir)/BCR -I$(srcdir)/BCRext -I$(srcdir)/countWords -I$(srcdir)/backtracker ${OPENMP_CXXFLAGS}

# Unfinished tool?
#indexBwt_SOURCES = BCRext/IndexBwt.cpp
#indexBwt_LDADD = liball.a
#indexBwt_CXXFLAGS = -I$(srcdir) -I$(srcdir)/shared -I$(srcdir)/BCR -I$(srcdir)/BCRext -I$(srcdir)/countWords


liball_a_SOURCES = \
	BCR/BCR_BWTCollection.cpp \
	BCR/BuildBCR.cpp \
	BCR/BWTCollection.cpp \
	BCR/BWTCollection.hh \
	BCR/BCRexternalBWT.cpp \
	BCR/BCRexternalBWT.hh \
	BCR/PredictiveEncoding.cpp \
	BCR/PredictiveEncoding.hh \
	BCR/Sorting.cpp \
	BCR/Sorting.hh \
	BCR/TransposeFasta.cpp \
	BCR/TransposeFasta.hh \
	shared/Algorithm.hh \
	shared/Alphabet.hh \
	shared/Config.hh \
	shared/EndPosFile.cpp \
	shared/EndPosFile.hh \
	shared/FragmentedVector.hh \
	shared/LetterCount.hh \
	shared/Types.hh \
	shared/Timer.cpp \
	shared/Timer.hh \
	shared/Tools.cpp \
	shared/Tools.hh \
	shared/Filename.cpp \
	shared/Filename.hh \
	shared/SeqReader.cpp \
	shared/SeqReader.hh \
	shared/SequenceExtractor.cpp \
	shared/SequenceExtractor.hh \
	parameters/BwtParameters.hh \
	parameters/CompareParameters.hh \
	parameters/ConvertParameters.hh \
	parameters/ExtendParameters.hh \
	parameters/IndexParameters.hh \
	parameters/SearchParameters.hh \
	parameters/UnbwtParameters.hh \
	BCRext/BwtReader.cpp \
	BCRext/BwtReader.hh \
	BCRext/BwtWriter.cpp \
	BCRext/BwtWriter.hh \
	BCRext/BwtIndex.cpp \
	BCRext/BwtIndex.hh \
	BCRext/ReadBuffer.cpp \
	BCRext/ReadBuffer.hh \
	BCRext/BCRext.cpp \
	BCRext/BCRext.hh \
	backtracker/BackTrackerBase.cpp \
	backtracker/BackTrackerBase.hh \
	backtracker/OneBwtBackTracker.cpp \
	backtracker/OneBwtBackTracker.hh \
	backtracker/TwoBwtBackTracker.cpp \
	backtracker/TwoBwtBackTracker.hh \
	backtracker/Range.cpp \
	backtracker/Range.hh \
	backtracker/RangeStore.cpp \
	backtracker/RangeStore.hh \
	backtracker/IntervalHandlerBase.cpp \
	backtracker/IntervalHandlerBase.hh \
	countWords/CountWords.cpp \
	countWords/CountWords.hh \
	countWords/IntervalHandlerMetagenome.cpp \
	countWords/IntervalHandlerMetagenome.hh \
	countWords/IntervalHandlerReference.cpp \
	countWords/IntervalHandlerReference.hh \
	countWords/IntervalHandlerSplice.cpp \
	countWords/IntervalHandlerSplice.hh \
	countWords/IntervalHandlerTumourNormal.cpp \
	countWords/IntervalHandlerTumourNormal.hh \
	errors/AlignmentParameters.hh \
	errors/BwtCorrectorParameters.hh \
	errors/WitnessReader.cpp \
	errors/WitnessReader.hh \
	errors/BwtCorrectorIntervalHandler.cpp \
	errors/BwtCorrectorIntervalHandler.hh \
	errors/BwtCorrector.cpp \
	errors/BwtCorrector.hh \
	errors/HiTECStats.cpp \
	errors/HiTECStats.hh \
	errors/CorrectionAligner.cpp \
	errors/CorrectionAligner.hh \
	errors/ErrorInfo.cpp \
	errors/ErrorInfo.hh \
	errors/ErrorCorrectionRange.cpp \
	errors/ErrorCorrectionRange.hh \
	search/SearchUsingBacktracker.cpp \
	search/SearchUsingBacktracker.hh \
	search/KmerSearchIntervalHandler.cpp \
	search/KmerSearchIntervalHandler.hh \
	search/KmerSearchRange.cpp \
	search/KmerSearchRange.hh \
	search/IntervalFile.cpp \
	search/IntervalFile.hh \
	search/Extender.cpp \
	search/Extender.hh \
	search/ExtenderIntervalHandler.cpp \
	search/ExtenderIntervalHandler.hh

liball_a_CXXFLAGS = -I$(srcdir) -I$(srcdir)/shared -I$(srcdir)/BCR -I$(srcdir)/BCRext -I$(srcdir)/countWords -I$(srcdir)/backtracker ${BOOST_CXXFLAGS} ${OPENMP_CXXFLAGS}


libzoo_a_SOURCES = \
	libzoo/cli/Common.cpp \
	libzoo/cli/Common.hh \
	libzoo/cli/ToolParameters.cpp \
	libzoo/cli/ToolParameters.hh \
	libzoo/io/Bcl.cpp \
	libzoo/io/Bcl.hh \
	libzoo/io/FastOFStream.cpp \
	libzoo/io/FastOFStream.hh \
	libzoo/util/Logger.cpp \
	libzoo/util/Logger.hh \
	libzoo/util/TemporaryFilesManager.cpp \
	libzoo/util/TemporaryFilesManager.hh \
	libzoo/util/AutoGrowVector.hh \
	libzoo/util/ColorText.cpp \
	libzoo/util/ColorText.hh

libzoo_a_CXXFLAGS = -I$(srcdir) ${BOOST_CXXFLAGS} ${OPENMP_CXXFLAGS}


# Unused file?
# Transpose.o: BCR/Transpose.cpp
