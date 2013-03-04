bin_PROGRAMS = beetl-bwt beetl-unbwt beetl-convert beetl-search beetl-compare
bin_SCRIPTS = beetl

beetl_bwt_SOURCES = BeetlBwt.cpp Common.cpp DatasetMetadata.cpp
beetl_bwt_LDADD = ../liball.a
beetl_bwt_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext

beetl_unbwt_SOURCES = BeetlUnbwt.cpp Common.cpp
beetl_unbwt_LDADD = ../liball.a
beetl_unbwt_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext

beetl_convert_SOURCES = BeetlConvert.cpp Common.cpp
beetl_convert_LDADD = ../liball.a
beetl_convert_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext

beetl_search_SOURCES = BeetlSearch.cpp Common.cpp
beetl_search_LDADD = ../liball.a
beetl_search_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext

beetl_compare_SOURCES = BeetlCompare.cpp Common.cpp
beetl_compare_LDADD = ../liball.a
beetl_compare_CXXFLAGS = -I$(srcdir)/.. -I$(srcdir)/../shared -I$(srcdir)/../BCR -I$(srcdir)/../BCRext
