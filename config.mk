# compiler and flags
CXX = g++
FLAGS = -g -Wno-deprecated -fPIC -fno-inline -Wno-write-strings
FLAGS += -fmax-errors=3
# extra flags
#FLAGS += -O0

# extra flags for Mac OS
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
	FLAGS += -std=c++11
endif

# ROOT
DEPS = -I$(shell root-config --incdir)
LIBS = $(shell root-config --glibs)
#LIBS += -lMinuit -lRooFitCore -lRooFit -lRooStats -lProof -lMathMore

# DELPHES
ifndef EXCLUDE_DELPHES
	DEPS += -I${DELPHES_HOME} -I${DELPHES_HOME}/external
	LIBS += -L${DELPHES_HOME} -lDelphes
	FLAGS += -DINCLUDE_DELPHES
endif

# DELPHES plugin: Fastjet Centauro
ifdef INCLUDE_CENTAURO
	LIBS += -L${DELPHES_HOME}/external/fastjet/plugins/Centauro -lCentauro
	DEPS += -I${DELPHES_HOME}/external/fastjet/plugins/Centauro
	FLAGS += -DINCLUDE_CENTAURO
endif

# MSTWPDF
DEPS += -I${MSTWPDF_HOME}
LIBS += -L${MSTWPDF_HOME} -lmstwpdf

# ADAGE
DEPS += -I${ADAGE_HOME}/include
LIBS += -L${ADAGE_HOME}/lib -lAdage

# SIDIS-EIC
ifndef SIDIS_EIC_HOME
	$(error "ERROR: run 'source environ.sh' before building")
endif
