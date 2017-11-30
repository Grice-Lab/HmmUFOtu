AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = -I$(top_srcdir)/include
noinst_PROGRAMS = bwt mksary sasearch suftest unbwt

bwt_SOURCES = bwt.c
bwt_LDADD = $(top_srcdir)/lib/libdivsufsort.a

mksary_SOURCES = mksary.c
mksary_LDADD = $(top_srcdir)/lib/libdivsufsort.a

sasearch_SOURCES = sasearch.c
sasearch_LDADD = $(top_srcdir)/lib/libdivsufsort.a

suftest_SOURCES = suftest.c
suftest_LDADD = $(top_srcdir)/lib/libdivsufsort.a

unbwt_SOURCES = unbwt.c
unbwt_LDADD = $(top_srcdir)/lib/libdivsufsort.a

