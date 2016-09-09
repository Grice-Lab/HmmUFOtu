SUBDIRS = libdivsufsort  libcds  math

AM_CXXFLAGS = -I$(top_srcdir)/src \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libcds/include

# Build a convenient library
noinst_LIBRARIES = libHmmUFOtu.a
libHmmUFOtu_a_SOURCES = \
EGMath.a \
HmmUFOtuConst.cpp \
CommandOptions.cpp \
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
CSFMIndex.cpp \
PhyloTree.cpp

# Build the executables
bin_PROGRAMS=hmmufotu hmmufotu-train-dm hmmufotu-train

hmmufotu_SOURCES=HmmUFOtu.cpp
hmmufotu_LDADD = libHmmUFOtu.a \
$(srcdir)/math/libEGMath.a
hmmufotu_LDFLAGS = -static

hmmufotu_train_dm_SOURCES=hmmufotu-train-dm.cpp
hmmufotu_train_dm_LDADD = libHmmUFOtu.a \
$(srcdir)/math/libEGMath.a
hmmufotu_train_dm_LDFLAGS = -static

hmmufotu_train_SOURCES = hmmufotu-train.cpp
hmmufotu_train_LDADD = libHmmUFOtu.a \
$(srcdir)/math/libEGMath.a
hmmufotu_train_CPPFLAGS = -DDATA=$(pkgdatadir)
hmmufotu_train_LDFLAGS = -static
