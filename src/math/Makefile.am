# This is the Makefile for the EGriceLab::Math library
AUTOMAKE_OPTIONS = subdir-objects

AM_CXXFLAGS = -I$(top_srcdir)/math

# Build a static/shared library
lib_LIBRARIES = libEGMath.a
libEGMath_a_SOURCES = Stats.cpp  LinearAlgebraBasic.cpp \
DirichletModel.cpp \
DirichletDensity.cpp \
DirichletMixture.cpp
