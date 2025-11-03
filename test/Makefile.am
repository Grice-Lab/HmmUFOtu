AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = \
-I$(top_srcdir) \
-I$(top_srcdir)/UCSC

AM_CXXFLAGS = \
-I$(top_srcdir) \
-I$(top_srcdir)/UCSC

check_PROGRAMS = \
StringUtils_test \
GFF_test

StringUtils_test_SOURCES = StringUtils_test.cpp
StringUtils_test_LDADD = \
$(top_srcdir)/libEGUtil.a

GFF_test_SOURCES = GFF_test.cpp
GFF_test_LDADD = \
$(top_srcdir)/UCSC/libEGUCSC.a \
$(top_srcdir)/libEGUtil.a

TESTS = \
StringUtils_test \
GFF_test
