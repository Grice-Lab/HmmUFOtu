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
AlphabetFactory.cpp \
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
TN93.cpp \
HKY85.cpp \
F81.cpp \
K80.cpp \
JC69.cpp \
DiscreteGammaModel.cpp \
DNASubModelFactory.cpp

#libHmmUFOtu_phylo_CPPFLAGS = $(AM_CFLAGS)
#libHmmUFOtu_phylo_CPPFLAGS += $(OPENMP_CFLAGS)
#libHmmUFOtu_phylo_CXXFLAGS = $(AM_CXXFLAGS)
#libHmmUFOtu_phylo_CXXFLAGS += $(OPENMP_CXXFLAGS)

# Build the executables
bin_PROGRAMS = \
hmmufotu \
hmmufotu-train-dm \
hmmufotu-train-sm \
hmmufotu-train-hmm \
hmmufotu-sim \
hmmufotu-place \
hmmufotu-build \
hmmufotu-inspect \
hmmufotu \
hmmufotu-sum

hmmufotu_train_dm_SOURCES=hmmufotu-train-dm.cpp
hmmufotu_train_dm_LDADD = libHmmUFOtu_common.a libHmmUFOtu_hmm.a math/libEGMath.a
hmmufotu_train_dm_LDFLAGS = -static

hmmufotu_train_sm_SOURCES=hmmufotu-train-sm.cpp
hmmufotu_train_sm_LDADD = libHmmUFOtu_phylo.a libHmmUFOtu_common.a math/libEGMath.a
hmmufotu_train_sm_LDFLAGS = -static

hmmufotu_train_hmm_SOURCES = hmmufotu-train-hmm.cpp
hmmufotu_train_hmm_LDADD = libHmmUFOtu_hmm.a math/libEGMath.a libHmmUFOtu_common.a
hmmufotu_train_hmm_LDFLAGS = -static
hmmufotu_train_hmm_CPPFLAGS = $(AM_CPPFLAGS) -DDM_DATADIR=\"$(pkgdatadir)\"

hmmufotu_sim_SOURCES = hmmufotu-sim.cpp
hmmufotu_sim_LDADD = libHmmUFOtu_phylo.a libHmmUFOtu_common.a
hmmufotu_sim_LDFLAGS = -static

hmmufotu_place_SOURCES = hmmufotu-place.cpp
hmmufotu_place_LDADD = libHmmUFOtu_phylo.a libHmmUFOtu_common.a
hmmufotu_place_LDFLAGS = -static

hmmufotu_build_SOURCES = hmmufotu-build.cpp
hmmufotu_build_LDADD = libHmmUFOtu_phylo.a libHmmUFOtu_hmm.a libHmmUFOtu_common.a math/libEGMath.a \
libdivsufsort/lib/libdivsufsort.a libcds/src/libcds.la
hmmufotu_build_LDFLAGS = -static
hmmufotu_build_CPPFLAGS = $(AM_CPPFLAGS) -DDM_DATADIR=\"$(pkgdatadir)\" -DSM_DATADIR=\"$(pkgdatadir)\"

hmmufotu_inspect_SOURCES = hmmufotu-inspect.cpp
hmmufotu_inspect_LDADD = libHmmUFOtu_phylo.a libHmmUFOtu_hmm.a libHmmUFOtu_common.a math/libEGMath.a \
libdivsufsort/lib/libdivsufsort.a libcds/src/libcds.la
hmmufotu_inspect_LDFLAGS = -static

hmmufotu_SOURCES = hmmufotu.cpp
hmmufotu_CPPFLAGS = $(AM_CFLAGS)
hmmufotu_CPPFLAGS += $(OPENMP_CFLAGS)
hmmufotu_CFLAGS = $(AM_CFLAGS)
hmmufotu_CFLAGS += $(OPENMP_CFLAGS)
hmmufotu_LDADD = libHmmUFOtu_phylo.a libHmmUFOtu_hmm.a math/libEGMath.a libHmmUFOtu_common.a \
libdivsufsort/lib/libdivsufsort.a libcds/src/libcds.la
hmmufotu_LDFLAGS = -static
hmmufotu_LDFLAGS += $(OPENMP_CFLAGS)

hmmufotu_sum_SOURCES = hmmufotu-sum.cpp
hmmufotu_sum_LDADD = libHmmUFOtu_phylo.a libHmmUFOtu_common.a math/libEGMath.a
hmmufotu_sum_LDFLAGS = -static
