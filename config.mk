# compiler and flags
CXX = g++
FLAGS = -Wno-deprecated -fPIC
FLAGS += -fmax-errors=3
# FLAGS += -fvisibility=hidden # FIXME: required by pybind, but causes unresolved symbols in cling...
# FLAGS += -g # build with debugging symbols
# FLAGS += -O0

# extra flags for Mac OS
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
	FLAGS += -std=c++11
endif

# ROOT
DEPS = -I$(shell root-config --incdir)
LIBS = $(shell root-config --glibs)
#LIBS += -lMinuit -lRooFitCore -lRooFit -lRooStats -lProof -lMathMore

# Data Model (PODIO + EDM4hep + EDM4eic)
LIBS += -L/usr/local/lib -lpodio -lpodioRootIO -ledm4hep -ledm4eic

# PYTHON (for pybind)
LIBS += -lpython3.10

# Miscellaneous
LIBS += -lfmt

# DELPHES
ifdef EXCLUDE_DELPHES
	FLAGS += -DEXCLUDE_DELPHES
else
	DEPS += -I${DELPHES_HOME} -I${DELPHES_HOME}/external
	LIBS += -L${DELPHES_HOME} -lDelphes
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

# EPIC-ANALYSIS
ifndef EPIC_ANALYSIS_HOME
$(error "ERROR: run 'source environ.sh' before building")
endif
