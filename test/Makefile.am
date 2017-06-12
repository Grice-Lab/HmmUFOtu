AUTOMAKE_OPTIONS = subdir-objects

AM_CXXFLAGS = -I$(top_srcdir)/src \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libcds/include

noinst_PROGRAMS = \
MSAIO_test \
bHmmPrior_IO_test \
bHmm_IO_test \
dna_model_IO_test \
FMIO_test \
PTU_IO_test

MSAIO_test_SOURCES = MSAIO_test.cpp
MSAIO_test_LDADD = $(top_srcdir)/src/libHmmUFOtu_common.a

bHmmPrior_IO_test_SOURCES = bHmmPrior_IO_test.cpp
bHmmPrior_IO_test_LDADD = $(top_srcdir)/src/libHmmUFOtu_hmm.a \
$(top_srcdir)/src/libHmmUFOtu_common.a $(top_srcdir)/src/math/libEGMath.a

bHmm_IO_test_SOURCES = bHmm_IO_test.cpp
bHmm_IO_test_LDADD = $(top_srcdir)/src/libHmmUFOtu_hmm.a \
$(top_srcdir)/src/libHmmUFOtu_common.a $(top_srcdir)/src/math/libEGMath.a

dna_model_IO_test_SOURCES = dna_model_IO_test.cpp
dna_model_IO_test_LDADD = $(top_srcdir)/src/libHmmUFOtu_phylo.a \
$(top_srcdir)/src/libHmmUFOtu_common.a $(top_srcdir)/src/math/libEGMath.a
dna_model_IO_test_LDFLAGS = -static

FMIO_test_SOURCES = FMIO_test.cpp
FMIO_test_LDADD = $(top_srcdir)/src/libHmmUFOtu_hmm.a \
$(top_srcdir)/src/libHmmUFOtu_common.a \
$(top_srcdir)/src/libdivsufsort/lib/libdivsufsort.a \
$(top_srcdir)/src/libcds/src/libcds.la
FMIO_test_LDFLAGS = -static

PTU_IO_test_SOURCES = PTU_IO_test.cpp
PTU_IO_test_LDADD = $(top_srcdir)/src/libHmmUFOtu_phylo.a \
$(top_srcdir)/src/libHmmUFOtu_common.a
PTU_IO_test_LDFLAGS = -static