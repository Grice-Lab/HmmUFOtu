SUBDIRS = libdivsufsort  libcds  math

AM_CXXFLAGS = -I$(top_srcdir)/src \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libcds/include

# Build convenient libraries
noinst_LIBRARIES = libHmmUFOtu_common.a libHmmUFOtu_hmm.a libHmmUFOtu_phylo.a

libHmmUFOtu_common_a_SOURCES = \
HmmUFOtuConst.cpp \
HmmUFOtuEnv.cpp \
ProgLog.cpp \
CommandOptions.cpp \
DegenAlphabet.cpp \
StringUtils.cpp \
IUPACNucl.cpp \
IUPACAmino.cpp \
DNA.cpp \
SeqCommons.cpp \
PrimarySeq.cpp \
DigitalSeq.cpp \
SeqIO.cpp \
MSA.cpp

libHmmUFOtu_hmm_a_SOURCES = \
BandedHMMP7Bg.cpp \
BandedHMMP7Prior.cpp \
BandedHMMP7.cpp \
CSFMIndex.cpp

libHmmUFOtu_phylo_a_SOURCES = \
NewickTree.cpp \
PhyloTreeUnrooted.cpp \
DNASubModel.cpp \
GTR.cpp \
DNASubModelFactory.cpp

# Build the executables
bin_PROGRAMS=hmmufotu hmmufotu-train-dm hmmufotu-train

hmmufotu_SOURCES=HmmUFOtu.cpp
hmmufotu_LDADD = libHmmUFOtu_common.a math/libEGMath.a
hmmufotu_LDFLAGS = -static

hmmufotu_train_dm_SOURCES=hmmufotu-train-dm.cpp
hmmufotu_train_dm_LDADD = libHmmUFOtu_common.a libHmmUFOtu_hmm.a math/libEGMath.a
hmmufotu_train_dm_LDFLAGS = -static

hmmufotu_train_SOURCES = hmmufotu-train.cpp
hmmufotu_train_LDADD = libHmmUFOtu_common.a libHmmUFOtu_hmm.a math/libEGMath.a
hmmufotu_train_CPPFLAGS = -DDATA=$(pkgdatadir)
hmmufotu_train_LDFLAGS = -static
