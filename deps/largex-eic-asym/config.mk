# compiler and flags
CXX = g++
FLAGS = -g -Wno-deprecated -fPIC -fno-inline -Wno-write-strings
FLAGS += -fmax-errors=3

# extra flags for Mac OS
#UNAME := $(shell uname)
#ifeq ($(UNAME), Darwin)
    #FLAGS += -std=c++11
#endif

# extra flags for valgrind
#FLAGS += -O0

# ROOT
DEPS = -I$(shell root-config --incdir)
LIBS = $(shell root-config --glibs)
LIBS += -lMinuit -lRooFitCore -lRooFit -lRooStats -lProof -lMathMore

# BruFit
DEPS += -I$(BRUFIT)/core
LIBS += -L$(BRUFIT)/lib -lbrufit

# shared object name and source directory
LARGEXASYM = LargexAsym
LARGEXASYMOBJ := lib$(LARGEXASYM).so
