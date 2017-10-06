SUBDIRS = libdivsufsort  libcds  math

AM_CFLAGS = -I$(top_srcdir)/src \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libcds/include \
$(BOOST_CPPFLAGS)

AM_CXXFLAGS = -I$(top_srcdir)/src \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libcds/include \
$(BOOST_CPPFLAGS)

AM_CPPFLAGS = $(BOOST_CPPFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS)

# Build convenient libraries
noinst_LIBRARIES = libHmmUFOtu_common.a libHmmUFOtu_hmm.a libHmmUFOtu_phylo.a libHmmUFOtu_OTU.a

libHmmUFOtu_common_a_SOURCES = \
ProgLog.cpp \
CommandOptions.cpp \
VersionSequence.cpp \
DegenAlphabet.cpp \
StringUtils.cpp \
IUPACNucl.cpp \
IUPACAmino.cpp \
DNA.cpp \
AlphabetFactory.cpp \
PrimarySeq.cpp \
DigitalSeq.cpp \
SeqIO.cpp \
SeqUtils.cpp \
MSA.cpp \
CSLoc.cpp \
TSVRecord.cpp \
TSVScanner.cpp

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

libHmmUFOtu_OTU_a_SOURCES = \
OTUObserved.cpp \
OTUTable.cpp

# Build the executables
bin_PROGRAMS = \
hmmufotu \
hmmufotu-train-dm \
hmmufotu-train-sm \
hmmufotu-train-hmm \
hmmufotu-sim \
hmmufotu-build \
hmmufotu-inspect \
hmmufotu \
hmmufotu-sum \
hmmufotu-anneal \
hmmufotu-subset \
hmmufotu-norm

hmmufotu_train_dm_SOURCES=hmmufotu-train-dm.cpp HmmUFOtuEnv.cpp
hmmufotu_train_dm_LDADD = libHmmUFOtu_common.a libHmmUFOtu_hmm.a math/libEGMath.a

hmmufotu_train_sm_SOURCES=hmmufotu-train-sm.cpp HmmUFOtuEnv.cpp
hmmufotu_train_sm_LDADD = libHmmUFOtu_phylo.a libHmmUFOtu_common.a math/libEGMath.a

hmmufotu_train_hmm_SOURCES = hmmufotu-train-hmm.cpp HmmUFOtuEnv.cpp
hmmufotu_train_hmm_LDADD = libHmmUFOtu_hmm.a math/libEGMath.a libHmmUFOtu_common.a
hmmufotu_train_hmm_CPPFLAGS = -DSRC_DATADIR=\"$(abs_top_srcdir)/data\" -DPKG_DATADIR=\"$(pkgdatadir)\"

hmmufotu_sim_SOURCES = hmmufotu-sim.cpp HmmUFOtuEnv.cpp
hmmufotu_sim_LDADD = libHmmUFOtu_phylo.a libHmmUFOtu_common.a

hmmufotu_build_SOURCES = hmmufotu-build.cpp HmmUFOtuEnv.cpp
hmmufotu_build_LDADD = libHmmUFOtu_phylo.a libHmmUFOtu_hmm.a libHmmUFOtu_common.a math/libEGMath.a \
libdivsufsort/lib/libdivsufsort.a libcds/src/libcds.la
hmmufotu_build_CPPFLAGS = -DSRC_DATADIR=\"$(abs_top_srcdir)/data\" -DPKG_DATADIR=\"$(pkgdatadir)\"

hmmufotu_inspect_SOURCES = hmmufotu-inspect.cpp HmmUFOtuEnv.cpp
hmmufotu_inspect_LDADD = libHmmUFOtu_phylo.a libHmmUFOtu_hmm.a libHmmUFOtu_common.a math/libEGMath.a \
libdivsufsort/lib/libdivsufsort.a libcds/src/libcds.la

hmmufotu_SOURCES = hmmufotu.cpp HmmUFOtu_main.cpp HmmUFOtuEnv.cpp
hmmufotu_CPPFLAGS = $(OPENMP_CFLAGS) $(BOOST_CPPFLAGS)
hmmufotu_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CFLAGS)
hmmufotu_LDADD = libHmmUFOtu_phylo.a libHmmUFOtu_hmm.a math/libEGMath.a libHmmUFOtu_common.a \
libdivsufsort/lib/libdivsufsort.a libcds/src/libcds.la \
$(BOOST_IOSTREAMS_LIB)
hmmufotu_LDFLAGS = $(BOOST_LDFLAGS)

hmmufotu_sum_SOURCES = hmmufotu-sum.cpp HmmUFOtu_main.cpp HmmUFOtuEnv.cpp
hmmufotu_sum_LDADD = libHmmUFOtu_OTU.a libHmmUFOtu_phylo.a libHmmUFOtu_hmm.a libHmmUFOtu_common.a math/libEGMath.a \
libdivsufsort/lib/libdivsufsort.a libcds/src/libcds.la

hmmufotu_anneal_SOURCES = hmmufotu-anneal.cpp HmmUFOtuEnv.cpp
hmmufotu_anneal_LDADD = libHmmUFOtu_phylo.a libHmmUFOtu_hmm.a libHmmUFOtu_common.a math/libEGMath.a \
libdivsufsort/lib/libdivsufsort.a libcds/src/libcds.la

hmmufotu_subset_SOURCES = hmmufotu-subset.cpp HmmUFOtuEnv.cpp
hmmufotu_subset_LDADD = libHmmUFOtu_OTU.a libHmmUFOtu_common.a

hmmufotu_norm_SOURCES = hmmufotu-norm.cpp HmmUFOtuEnv.cpp
hmmufotu_norm_LDADD = libHmmUFOtu_OTU.a libHmmUFOtu_common.a
