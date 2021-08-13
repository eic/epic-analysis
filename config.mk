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
DEPS += -I${DELPHES_HOME} -I${DELPHES_HOME}/external
LIBS += -L${DELPHES_HOME} -lDelphes 

INCCENTAURO = 
ifdef INCCENTAURO
LIBS+= -L${DELPHES_HOME}/external/fastjet/plugins/Centauro -lCentauro
DEPS+= -I${DELPHES_HOME}/external/fastjet/plugins/Centauro
endif

# shared object name and source directory
LARGEX = Largex
LARGEXOBJ := lib$(LARGEX).so
