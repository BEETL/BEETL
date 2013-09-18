AM_CXXFLAGS = ${OPENMP_CXXFLAGS}
AM_LDFLAGS = -L${BOOST_ROOT}/lib

bin_PROGRAMS = beetl-bwt beetl-unbwt beetl-convert beetl-search beetl-compare beetl-correct align-corrector-strings beetl-hitec
bin_SCRIPTS = beetl

beetl_bwt_SOURCES = BeetlBwt.cpp Common.cpp DatasetMetadata.cpp
beetl_bwt_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD}
beetl_bwt_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext ${OPENMP_CXXFLAGS}

beetl_unbwt_SOURCES = BeetlUnbwt.cpp
beetl_unbwt_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD}
beetl_unbwt_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext ${OPENMP_CXXFLAGS}

beetl_convert_SOURCES = BeetlConvert.cpp
beetl_convert_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD} -lz
beetl_convert_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext ${OPENMP_CXXFLAGS}

beetl_search_SOURCES = BeetlSearch.cpp
beetl_search_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD}
beetl_search_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext ${OPENMP_CXXFLAGS}

beetl_compare_SOURCES = BeetlCompare.cpp Common.cpp
beetl_compare_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD}
beetl_compare_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext ${OPENMP_CXXFLAGS}

beetl_correct_SOURCES = BeetlCorrect.cpp
beetl_correct_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD}
beetl_correct_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext ${OPENMP_CXXFLAGS}

align_corrector_strings_SOURCES = AlignCorrectorStrings.cpp
align_corrector_strings_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD}
align_corrector_strings_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext ${OPENMP_CXXFLAGS}

beetl_hitec_SOURCES = BeetlHiTEC.cpp
beetl_hitec_LDADD = ../liball.a ../libzoo.a ${BOOST_LDADD}
beetl_hitec_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext ${OPENMP_CXXFLAGS}
