# This is the Makefile for the EGriceLab::Math library
ACLOCAL_AMFLAGS = -I m4

AM_CXXFLAGS = $(BOOST_CPPFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS)

# Build a static/shared library
lib_LIBRARIES = libEGMath.a
libEGMath_a_SOURCES = \
Stats.cpp \
LinearAlgebraBasic.cpp \
DirichletModel.cpp \
DirichletDensity.cpp \
DirichletMixture.cpp \
RootFinder.cpp