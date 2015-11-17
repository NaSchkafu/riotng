ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = @IV_CPPFLAGS@ @DEBUG_FLAGS@
lib_LTLIBRARIES = libtm.la
libtm_la_SOURCES = tm-base/taylormodel.cpp tm-base/adaptintval.cpp tm-base/polyeval.cpp tm-base/taylormodel_impl.cpp tm-base/primetable.cpp \
	tm-base/monom.cpp tm-base/factorialtable.cpp
nobase_include_HEADERS = tm-base/actions.h tm-base/factorialtable.h tm-base/observer.h tm-base/hashfunc.h tm-base/taylormodel.h \
	tm-base/adaptintval.h tm-base/hashtable.h tm-base/polyeval.h tm-base/degreeeval.h tm-base/taylormodel_impl.h \
	tm-base/error.h tm-base/primetable.h tm-base/vardatatable.h tm-base/monom.h tm-base/refcounter.h tm-base/vartable.h
libtm_la_LIBADD = @IV_LIBS@
libtm_la_LDFLAGS = -version-info 0:1:0

if ODE_SUPPORT
bin_PROGRAMS = riotng
riotng_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/src/tm-base
riotng_SOURCES = dgl/dvector.cpp dgl/linsys.cpp dgl/matrix.cpp dgl/data.cpp dgl/imatrix.cpp dgl/main.cpp \
	dgl/vvtaylormodel.cpp dgl/shrinkwrap.cpp dgl/dglalg.cpp dgl/ivector.cpp dgl/matinv.cpp \
	dgl/stopwatch.h dgl/dvector.h dgl/linsys.h dgl/matrix.h dgl/timeid.h \
	dgl/data.h dgl/output.h dgl/degree.h dgl/imatrix.h dgl/vvtaylormodel.h \
	dgl/shrinkwrap.h dgl/dglalg.h dgl/ivector.h dgl/matinv.h

riotng_LDADD = libtm.la
endif

