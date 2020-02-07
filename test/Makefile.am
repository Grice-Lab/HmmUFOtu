AUTOMAKE_OPTIONS = subdir-objects

AM_CXXFLAGS = -I$(top_srcdir)/src \
-I$(top_srcdir)/src/util \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libcds/include \
$(BOOST_CPPFLAGS) \
$(OPENMP_CFLAGS)

AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(OPENMP_CFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS) $(OPENMP_CFLAGS)

check_PROGRAMS = \
MSAIO_test \
bHmmPrior_IO_test \
bHmm_IO_test \
dna_model_IO_test \
FMIO_test \
PTU_IO_test \
CSFMIndex_test

MSAIO_test_SOURCES = MSAIO_test.cpp
MSAIO_test_LDADD = $(top_srcdir)/src/libHmmUFOtu_common.a $(top_srcdir)/src/util/libEGUtil.a \
$(top_srcdir)/src/HmmUFOtuEnv.o

bHmmPrior_IO_test_SOURCES = bHmmPrior_IO_test.cpp
bHmmPrior_IO_test_LDADD = $(top_srcdir)/src/libHmmUFOtu_hmm.a $(top_srcdir)/src/libHmmUFOtu_common.a \
$(top_srcdir)/src/util/libEGUtil.a $(top_srcdir)/src/math/libEGMath.a \
$(top_srcdir)/src/HmmUFOtuEnv.o

bHmm_IO_test_SOURCES = bHmm_IO_test.cpp
bHmm_IO_test_LDADD = $(top_srcdir)/src/libHmmUFOtu_hmm.a $(top_srcdir)/src/libHmmUFOtu_common.a \
$(top_srcdir)/src/util/libEGUtil.a $(top_srcdir)/src/math/libEGMath.a \
$(top_srcdir)/src/HmmUFOtuEnv.o

dna_model_IO_test_SOURCES = dna_model_IO_test.cpp
dna_model_IO_test_LDADD = $(top_srcdir)/src/libHmmUFOtu_phylo.a $(top_srcdir)/src/libHmmUFOtu_common.a \
$(top_srcdir)/src/util/libEGUtil.a $(top_srcdir)/src/math/libEGMath.a \
$(top_srcdir)/src/HmmUFOtuEnv.o

FMIO_test_SOURCES = FMIO_test.cpp
FMIO_test_LDADD = $(top_srcdir)/src/libHmmUFOtu_hmm.a $(top_srcdir)/src/libHmmUFOtu_common.a \
$(top_srcdir)/src/util/libEGUtil.a \
$(top_srcdir)/src/libdivsufsort/lib/libdivsufsort.a \
$(top_srcdir)/src/libcds/src/libcds.la \
$(top_srcdir)/src/HmmUFOtuEnv.o

PTU_IO_test_SOURCES = PTU_IO_test.cpp
PTU_IO_test_LDADD = $(top_srcdir)/src/libHmmUFOtu_phylo.a $(top_srcdir)/src/libHmmUFOtu_common.a \
$(top_srcdir)/src/util/libEGUtil.a \
$(top_srcdir)/src/math/libEGMath.a \
$(top_srcdir)/src/HmmUFOtuEnv.o

CSFMIndex_test_SOURCES = CSFMIndex_test.cpp
CSFMIndex_test_LDADD = $(top_srcdir)/src/libHmmUFOtu_hmm.a $(top_srcdir)/src/libHmmUFOtu_common.a \
$(top_srcdir)/src/util/libEGUtil.a \
$(top_srcdir)/src/math/libEGMath.a \
$(top_srcdir)/src/libdivsufsort/lib/libdivsufsort.a \
$(top_srcdir)/src/libcds/src/libcds.la \
$(top_srcdir)/src/HmmUFOtuEnv.o

TESTS = CSFMIndex_test GTR-t.sh TN93-t.sh HKY85-t.sh GTR-dG-t.sh 
if HAVE_JSONCPP
TESTS += jplace-t.sh
endif
TESTS += sim-run-SE-t.sh 
TESTS += sim-run-PE-t.sh
