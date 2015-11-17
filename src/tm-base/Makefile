# Which C++ compiler? 
# Possible values: gnu         for the GNU C/C++ Compiler (gcc)
#                  intel       for the INTEL C/C++ Compiler (OS Linux)
COMPILER=gnu

# Which library and path (compiled with the above compiler!).
# Possible values: filib++ 
#                  cxsc    
LIBRARY=cxsc# NOT USED YET!!!!!
DIR=/usr/local

ifeq ($(COMPILER),gnu)
# settings for GNU - Compiler ----------------------------------------
CC=gcc#                  name of the C compiler
CCOPTS=-Wall -frounding-math -ffloat-store#            optional flags forwarded to the C compiler
CCOPTIMIZE=#             -O0/-O1/-O2/-O3  flag for optimization
CXX=g++#                 name of the C++ compiler
CXXOPTS=#                optional flags forwarded to the C++ compiler
CXXOPTIMIZE=-O3#         -O0/-O1/-O2/-O3 flag for optimization
else
ifeq ($(COMPILER),intel)
# settings for Intel C++ Compiler ------------------------------------
CC=icc#                  name of the C compiler
CCOPTS=#                 optional flags forwarded to the C compiler
CCOPTIMIZE=#             -O0/-O1/-O2/-O3  flag for optimization
CXX=icpc#                name of the C++ compiler
CXXOPTS=-w#              optional flags forwarded to the C++ compiler
CXXOPTIMIZE=-O3#         -O0/-O1/-O2/-O3 flag for optimization
endif
endif

# Settings for all systems -------------------------------------------
#
SOURCES	= adaptintval.cpp polyeval.cpp taylormodel_impl.cpp \
          factorialtable.cpp monom.cpp primetable.cpp taylormodel.cpp
#
AR = ar
INCLUDE = -I$(DIR)/include/
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
	echo "Zur Verfügung stehende Makeziele:"
	echo
	echo "lib-filib:  Erzeugen der TM-Bibliothek basierend auf filib++."
	echo "lib-cxsc :  Erzeugen der TM-Bibliothek basierend auf cxsc."
	echo "clean:      Löschen aller temporären Dateien."
	echo
	echo "--------------------------------------------------------------------------------"

# Creating the Library based on filib++.
PHONY:   lib-filib
.SILENT: lib-filib
lib-filib: CXXOPTS += -DFILIB_VERSION -frounding-math -ffloat-store
lib-filib: TEXT = filib++
lib-filib: clean print_first lib print_last
	mv libtm.a libtmfilib.a

# Creating the Library based on C-XSC.
PHONY:   lib-cxsc
.SILENT: lib-cxsc
lib-cxsc: CXXOPTS += -DCXSC_VERSION
lib-cxsc: TEXT = C-XSC
lib-cxsc: clean print_first lib print_last
	mv libtm.a libtmcxsc.a

# ------------------------------------------------------------------------------
# Creating the library
#
PHONY: lib

$(OBJECTS): %.o: %.cpp
	$(CXX) $(CXXOPTIMIZE) $(CXXOPTS) $(INCLUDE) -c $< -o $@

lib: libtm.a

libtm.a: $(OBJECTS)
	$(AR) rs $@ $?

# Header Files
ACTIONS.H = actions.h $(OBSERVER.H)
ADAPTINTVAL.H = adaptintval.h
DEGREEEVAL.H = degreeeval.h $(OBSERVER.H) $(ACTIONS.H) $(VARTABLE.H) $(MONOM.H) 
ERROR.H = error.h
FACTORIALTABLE.H = factorialtable.h $(ADAPTINTVAL.H)
HASHFUNC.H = hashfunc.h $(HASHTABLE.H)
HASHTABLE.H = hashtable.h $(PRIMETABLE.H)
MONOM.H = monom.h $(VARTABLE.H)
OBSERVER.H = observer.h
POLYEVAL.H = polyeval.h $(ADAPTINTVAL.H) $(HASHFUNC.H) $(REFCOUNTER.H) $(VARDATATABLE.H)
PRIMETABLE.H = primetable.h
REFCOUNTER.H = refcounter.h
TAYLORMODEL.H = taylormodel.h $(REFCOUNTER.H) $(TAYLORMODEL_IMPL.H)
TAYLORMODEL_IMPL.H = taylormodel_impl.h $(HASHFUNC.H) $(ADAPTINTVAL.H) $(MONOM.H) $(VARDATATABLE.H) $(POLYEVAL.H) $(REFCOUNTER.H)
VARDATATABLE.H = vardatatable.h $(ADAPTINTVAL.H)
VARTABLE.H = vartable.h $(OBSERVER.H) $(ACTIONS.H)

# Object-Files
adaptintval.o: adaptintval.cpp
taylormodel.o: taylormodel.cpp $(TAYLORMODEL.H)
monom.o: monom.cpp $(MONOM.H)
polyeval.o: polyeval.cpp $(POLYEVAL.H)
taylormodel_impl.o: taylormodel_impl.cpp $(TAYLORMODEL_IMPL.H) $(QUICKSORT.H)
primetable.o: primetable.cpp $(PRIMETABLE.H)
factorialtable.o: factorialtable.cpp $(FACTORIALTABLE.H)

#
# ------------------------------------------------------------------------------
.PHONY:  print_first
.SILENT: print_first
print_first:
	echo
	echo "Erzeugen der TM-Bibliothek basierend auf" $(TEXT)
	echo
	echo 
	echo

.PHONY:  print_last
.SILENT: print_last
print_last:
	echo
	echo "Fertig!"
	echo


#
# Cleaning the directory
.PHONY:  clean
.SILENT: clean
clean: 
	find . -name "*~"    -print0 | xargs -0 rm 2>/dev/null; echo Removing temporary files; true
	find . -name "*.o"   -print0 | xargs -0 rm 2>/dev/null; echo Removing object files; true
