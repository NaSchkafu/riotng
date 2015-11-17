# C++ compiler.
CXX = g++
# Path to interval library.
DIR = /usr/local
#
#-------------------------------------------------------------------------------
#
TMDIR = ../tm-base
#
SOURCES	= matinv.cpp linsys.cpp matrix.cpp imatrix.cpp vvtaylormodel.cpp \
	ivector.cpp dvector.cpp shrinkwrap.cpp dglalg.cpp data.cpp main.cpp
#
OPTIONS = -O3
INCLUDE = -I$(DIR)/include -I$(TMDIR)
LIB = -L$(DIR)/lib -L$(TMDIR)
LINKERFLAGS = $(LIBTMFLAG)
#
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))
#
#-------------------------------------------------------------------------------
#
start: info

# Information
PHONY:   info
.SILENT: info
info:
	echo "--------------------------------------------------------------------------------"
	echo
	echo "Make targets:"
	echo
	echo "solve-filib:  Compile RiOT using filib++."
	echo "solve-cxsc :  Compile RiOT using C-XSC."
	echo "clean      :  Delete all temporary Files."
	echo
	echo "--------------------------------------------------------------------------------"

# Creating RiOT based on filib++.
PHONY:   solve-filib
solve-filib: OPTIONS    += -DFILIB_VERSION -frounding-math -ffloat-store
solve-filib: LIB        += -Wl,-rpath=$(DIR)/lib
solve-filib: LINKERFLAGS = -static -ltmfilib -lprim -lm 
solve-filib: $(TMDIR)/libtmfilib.a solve
	mv solve solve_filib

$(TMDIR)/libtmfilib.a: $(TMDIR)/Makefile $(TMDIR)/*.cpp $(TMDIR)/*.h
	@make -C $(TMDIR) DIR=$(DIR) lib-filib;

# Creating RiOT based on C-XSC.
PHONY:   solve-cxsc
solve-cxsc: OPTIONS    += -DCXSC_VERSION
solve-cxsc: LINKERFLAGS = -static -ltmcxsc -lcxsc 
solve-cxsc: $(TMDIR)/libtmcxsc.a solve
	mv solve solve_cxsc

$(TMDIR)/libtmcxsc.a: $(TMDIR)/Makefile $(TMDIR)/*.cpp $(TMDIR)/*.h
	@make -C $(TMDIR) DIR=$(DIR) lib-cxsc;

# ------------------------------------------------------------------------------
# Creating RiOT
#
PHONY: solve
solve: $(OBJECTS)
	$(CXX) $(LIB) $(OBJECTS) $(LINKERFLAGS) -o $@

$(OBJECTS): %.o: %.cpp
	$(CXX) $(OPTIONS) $(INCLUDE) -c $< -o $@

# Header Files
DATA.H = data.h
DEGREE.H = degree.h
DGLALG.H = dglalg.h $(DEGREE.H)
IMATRIX.H = imatrix.h
IVECTOR.H = ivector.h
LINSYS.H = linsys.h
MATRIX.H = matrix.h
MATINV.H = matinv.h
OUTPUT.H = output.h
SHRINKWRAP.H = shrinkwrap.h
TIMEID.H = timeid.h
DVECTOR.H = dvector.h
VVTAYLORMODEL.H = vvtaylormodel.h

# Object-Files
dglalg.o: dglalg.cpp $(DGLALG.H) $(SHRINKWRAP.H) $(IVECTOR.H) $(VVTAYLORMODEL.H) $(TIMEID.H) $(OUTPUT.H)
data.o: data.cpp  $(VVTAYLORMODEL.H) $(IVECTOR.H)
ivector.o: ivector.cpp $(IVECTOR.H)
main.o: main.cpp $(DGLALG.H) $(DATA.H) $(LIBDEP)
matinv.o: matinv.cpp $(MATINV.H) $(MATRIX.H)
matrix.o: matrix.cpp $(MATRIX.H)
imatrix.o: imatrix.cpp $(IMATRIX.H) $(LINSYS.H)
linsys.o: linsys.cpp $(LINSYS.H) $(MATINV.H) $(MATRIX.H) $(IVECTOR.H)
shrinkwrap.o: shrinkwrap.cpp $(SHRINKWRAP.H) $(IVECTOR.H) $(TIMEID.H)
dvector.o: dvector.cpp $(VECTOR.H)
vvtaylormodel.o: vvtaylormodel.cpp $(VVTAYLORMODEL.H) $(MATRIX.H) $(DVECTOR.H) $(IVECTOR.H)

#
# ------------------------------------------------------------------------------
# Cleaning the directory
.PHONY:  clean
.SILENT: clean
clean:
	find . -name "*~"    -print0 | xargs -0 rm 2>/dev/null; echo Removing temporary files; true
	find . -name "*.o"   -print0 | xargs -0 rm 2>/dev/null; echo Removing object files; true
