# compiler and flags
CXX = g++
FLAGS = -g -Wno-deprecated -fPIC -fno-inline -Wno-write-strings
FLAGS += -fmax-errors=3
#-fvisibility=hidden
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

# Data Model (PODIO + EDM4hep + EDM4eic)
LIBS += -L/usr/local/lib -lpodio -lpodioRootIO -ledm4hep -ledm4eic
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

#PYBIND
DEPS += -I/usr/._local/4lciz2i3o3tp525sudm44xu6vhqszeil/include/python3.10 -I/usr/._local/4lciz2i3o3tp525sudm44xu6vhqszeil/include/python3.10
LIBS += -L/usr/._local/4lciz2i3o3tp525sudm44xu6vhqszeil/lib -lcrypto -lintl -ldl  -lm -lm -lpython3.10
DEPS += -I/w/hallb-scshelf2102/clas12/pecar/eic/sidis-eic/deps/pybind11/include

# SIDIS-EIC
ifndef SIDIS_EIC_HOME
$(error "ERROR: run 'source environ.sh' before building")
endif
