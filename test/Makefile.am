AUTOMAKE_OPTIONS = subdir-objects

AM_CXXFLAGS = -I$(top_srcdir)/src

noinst_PROGRAMS = alphabet_test1 dna_test1   nucl_test1 amino_test1   \
seqIO_test1    banded_align_known_seed_test1 \
MSAIO_test1    MSAIO_test2

alphabet_test1_SOURCES = alphabet_test1.cpp \
../src/Alphabet.cpp \
../src/StringUtils.cpp

dna_test1_SOURCES = dna_test1.cpp \
../src/Alphabet.cpp \
../src/StringUtils.cpp \
../src/DegenAlphabet.cpp \
../src/DNA.cpp

nucl_test1_SOURCES = nucl_test1.cpp \
../src/Alphabet.cpp \
../src/StringUtils.cpp \
../src/DegenAlphabet.cpp \
../src/IUPACNucl.cpp

amino_test1_SOURCES = amino_test1.cpp \
../src/Alphabet.cpp \
../src/StringUtils.cpp \
../src/DegenAlphabet.cpp \
../src/IUPACAmino.cpp

seqIO_test1_SOURCES = seqIO_test1.cpp \
../src/Alphabet.cpp \
../src/StringUtils.cpp \
../src/DegenAlphabet.cpp \
../src/IUPACNucl.cpp \
../src/IUPACAmino.cpp \
../src/SeqCommons.cpp \
../src/PrimarySeq.cpp \
../src/SeqIO.cpp

banded_align_known_seed_test1_SOURCES = banded_align_known_seed_test1.cpp \
../src/Alphabet.cpp \
../src/StringUtils.cpp \
../src/DegenAlphabet.cpp \
../src/IUPACNucl.cpp \
../src/IUPACAmino.cpp \
../src/SeqCommons.cpp \
../src/PrimarySeq.cpp \
../src/BandedHMMP7Bg.cpp \
../src/BandedHMMP7.cpp

MSAIO_test1_SOURCES = MSAIO_test1.cpp \
../src/Alphabet.cpp \
../src/StringUtils.cpp \
../src/DegenAlphabet.cpp \
../src/IUPACNucl.cpp \
../src/IUPACAmino.cpp \
../src/SeqCommons.cpp \
../src/MSA.cpp \
../src/PrimarySeq.cpp \
../src/SeqIO.cpp \
../src/Stats.cpp

MSAIO_test2_SOURCES = MSAIO_test2.cpp \
../src/Alphabet.cpp \
../src/StringUtils.cpp \
../src/DegenAlphabet.cpp \
../src/IUPACNucl.cpp \
../src/IUPACAmino.cpp \
../src/SeqCommons.cpp \
../src/MSA.cpp \
../src/PrimarySeq.cpp \
../src/SeqIO.cpp \
../src/Stats.cpp
