SUBDIRS = libdivsufsort  libcds  math

AM_CXXFLAGS = -I$(top_srcdir)/src \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libcds/include

# Build a convenient library
noinst_LIBRARIES = libHmmUFOtu.a
libHmmUFOtu_a_SOURCES = \
EGMath.a \
DegenAlphabet.cpp \
StringUtils.cpp \
IUPACNucl.cpp \
IUPACAmino.cpp \
DNA.cpp \
SeqCommons.cpp \
PrimarySeq.cpp \
SeqIO.cpp \
BandedHMMP7Bg.cpp \
BandedHMMP7.cpp \
MSA.cpp \
CSFMIndex.cpp

# Build the executables
bin_PROGRAMS=HmmUFOtu
HmmUFOtu_SOURCES=HmmUFOtu.cpp HmmUFOtu.a

