AUTOMAKE_OPTIONS = subdir-objects

AM_CXXFLAGS = -I$(top_srcdir)/src \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libcds/include

noinst_PROGRAMS = alphabet_test1 dna_test1   nucl_test1 amino_test1   \
seqIO_test1    banded_align_known_seed_test1 \
MSAIO_test1    MSAIO_test2 \
FMIO_test1     FMIO_test2

alphabet_test1_SOURCES = alphabet_test1.cpp
alphabet_test1_LDADD = $(top_srcdir)/src/libHmmUFOtu.a

dna_test1_SOURCES = dna_test1.cpp
dna_test1_LDADD = $(top_srcdir)/src/libHmmUFOtu.a

nucl_test1_SOURCES = nucl_test1.cpp
nucl_test1_LDADD = $(top_srcdir)/src/libHmmUFOtu.a

amino_test1_SOURCES = amino_test1.cpp
amino_test1_LDADD = $(top_srcdir)/src/libHmmUFOtu.a

seqIO_test1_SOURCES = seqIO_test1.cpp
seqIO_test1_LDADD = $(top_srcdir)/src/libHmmUFOtu.a

banded_align_known_seed_test1_SOURCES = banded_align_known_seed_test1.cpp
banded_align_known_seed_test1_LDADD = $(top_srcdir)/src/libHmmUFOtu.a

MSAIO_test1_SOURCES = MSAIO_test1.cpp
MSAIO_test1_LDADD = $(top_srcdir)/src/libHmmUFOtu.a

MSAIO_test2_SOURCES = MSAIO_test2.cpp
MSAIO_test2_LDADD = $(top_srcdir)/src/libHmmUFOtu.a

FMIO_test1_SOURCES = FMIO_test1.cpp
FMIO_test1_LDADD = $(top_srcdir)/src/libHmmUFOtu.a \
$(top_srcdir)/src/libdivsufsort/lib/libdivsufsort.a \
$(top_srcdir)/src/libcds/src/libcds.la
FMIO_test1_LDFLAGS = -static

FMIO_test2_SOURCES = FMIO_test2.cpp
FMIO_test2_LDADD = $(top_srcdir)/src/libHmmUFOtu.a \
$(top_srcdir)/src/libdivsufsort/lib/libdivsufsort.a \
$(top_srcdir)/src/libcds/src/libcds.la
FMIO_test2_LDFLAGS = -static