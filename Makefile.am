SUBDIRS=. UCSC test

# This is the Makefile for the EGUtil library
ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = $(BOOST_CPPFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS)

# Build a static/shared library
lib_LIBRARIES = libEGUtil.a
libEGUtil_a_SOURCES = \
CommandOptions.cpp \
ProgLog.cpp \
StringUtils.cpp \
TSVRecord.cpp \
TSVScanner.cpp \
VersionSequence.cpp \
ProgEnv.cpp
