# Makefile for MSTW Fortran and C++ code.
# Comments to Graeme.Watt(at)cern.ch.

# If using g77, uncomment these two lines:
#FC=g77
#LDLIBS=-lg2c

# If using gfortran, uncomment these two lines:
FC=gfortran
LDLIBS=-lgfortran

# C and C++ compilers:
CC=gcc
CXX=g++

# Compiler flags:
FFLAGS=-Wall -W -fPIC
CXXFLAGS=-Wall -W -fPIC

all: libmstwpdf.so
debug: CXXFLAGS += -g
debug: clean all
release: CXXFLAGS += -O3
release: clean all

libmstwpdf.so : mstwpdf_cpp.o alphaS.o
	$(CXX) -shared $(CXXFLAGS) $< -o $@

mstwpdf.o : mstwpdf.f
	$(FC) -c $(FFLAGS) $< -o $@

alphaS.o : alphaS.f
	$(FC) -Wno-compare-reals -c $(FFLAGS) $< -o $@

#####################################
# C++ version:

mstwpdf_cpp.o : mstwpdf.cc mstwpdf.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

clean :
	-rm -f *.o libmstwpdf.so
