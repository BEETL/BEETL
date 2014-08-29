AM_CXXFLAGS = ${OPENMP_CXXFLAGS}
AM_LDFLAGS = -L${BOOST_ROOT}/lib

bin_PROGRAMS = beetl-bwt beetl-unbwt beetl-convert beetl-search beetl-compare beetl-correct beetl-correct-apply-corrections beetl-index beetl-extend
bin_SCRIPTS = beetl

beetl_bwt_SOURCES = BeetlBwt.cpp BeetlBwt.hh Common.cpp DatasetMetadata.cpp
beetl_bwt_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD}
beetl_bwt_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext ${OPENMP_CXXFLAGS}

beetl_unbwt_SOURCES = BeetlUnbwt.cpp BeetlUnbwt.hh
beetl_unbwt_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD}
beetl_unbwt_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext ${OPENMP_CXXFLAGS}

beetl_convert_SOURCES = BeetlConvert.cpp BeetlConvert.hh
beetl_convert_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD}
beetl_convert_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext ${OPENMP_CXXFLAGS}

beetl_search_SOURCES = BeetlSearch.cpp BeetlSearch.hh
beetl_search_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD}
beetl_search_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext ${OPENMP_CXXFLAGS}

beetl_compare_SOURCES = BeetlCompare.cpp BeetlCompare.hh Common.cpp
beetl_compare_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD}
beetl_compare_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext -I$(srcdir)/../backtracker ${OPENMP_CXXFLAGS}

beetl_correct_SOURCES = BeetlCorrect.cpp BeetlCorrect.hh
beetl_correct_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD}
beetl_correct_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext -I$(srcdir)/../backtracker ${OPENMP_CXXFLAGS}

beetl_correct_apply_corrections_SOURCES = AlignCorrectorStrings.cpp AlignCorrectorStrings.hh
beetl_correct_apply_corrections_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD}
beetl_correct_apply_corrections_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext -I$(srcdir)/../backtracker ${OPENMP_CXXFLAGS}

beetl_index_SOURCES = BeetlIndex.cpp BeetlIndex.hh
beetl_index_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD}
beetl_index_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext ${OPENMP_CXXFLAGS}

beetl_extend_SOURCES = BeetlExtend.cpp BeetlExtend.hh
beetl_extend_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD}
beetl_extend_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext -I$(srcdir)/../backtracker ${OPENMP_CXXFLAGS}


noinst_HEADERS = Common.hh DatasetMetadata.hh
