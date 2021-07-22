# compiler and flags
CXX = g++
FLAGS = -g -Wno-deprecated -fPIC -fno-inline -Wno-write-strings

# extra flags for Mac OS
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    FLAGS += -std=c++11
endif

# extra flags for valgrind
#FLAGS += -O0

# ROOT
DEPS = -I$(shell root-config --incdir)
LIBS = $(shell root-config --glibs)
#LIBS += -lMinuit -lRooFitCore -lRooFit -lRooStats -lProof -lMathMore

# DELPHES
DEPS += -I${DELPHES_HOME}
LIBS += -L${DELPHES_HOME} -lDelphes

# shared object name and source directory
LARGEX = Largex
LARGEXOBJ := lib$(LARGEX).so
