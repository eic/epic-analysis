# project name
PROJECT = EpicAnalysis

# check environment
ifndef EPIC_ANALYSIS_HOME
$(error "ERROR: run 'source environ.sh' before building")
endif


# compiler and flags
#################################################################

CXX := g++
FLAGS = -Wno-deprecated -fPIC
FLAGS += -fmax-errors=3
# FLAGS += -fvisibility=hidden # FIXME: required by pybind, but causes unresolved symbols in cling...
ROOTCLING = rootcling


# dependencies
#################################################################

# ROOT
DEP_INCLUDES  = -I$(shell root-config --incdir)
DEP_LIBRARIES = $(shell root-config --glibs)
#DEP_LIBRARIES += -lMinuit -lRooFitCore -lRooFit -lRooStats -lProof -lMathMore

# Data Model (PODIO + EDM4hep + EDM4eic)
DEP_LIBRARIES += -L/usr/local/lib -ledm4hep -ledm4eic -lpodio -lpodioRootIO

# Miscellaneous
DEP_LIBRARIES += -lfmt

# Delphes
ifdef EXCLUDE_DELPHES
	# optionally build without Delphes (for CI speedup)
	FLAGS += -DEXCLUDE_DELPHES
else
	DEP_INCLUDES  += -I${DELPHES_HOME} -I${DELPHES_HOME}/external
	DEP_LIBRARIES += -L${DELPHES_HOME} -lDelphes
endif

# Delphes plugin: Fastjet Centauro
ifdef INCLUDE_CENTAURO
	DEP_INCLUDES  += -I${DELPHES_HOME}/external/fastjet/plugins/Centauro
	DEP_LIBRARIES += -L${DELPHES_HOME}/external/fastjet/plugins/Centauro -lCentauro
	FLAGS += -DINCLUDE_CENTAURO
endif

# MSTWPDF
DEP_INCLUDES  += -I${MSTWPDF_HOME}
DEP_LIBRARIES += -L${MSTWPDF_HOME} -lmstwpdf
DEP_TARGETS   =  ${MSTWPDF_HOME}/libmstwpdf.so

# ADAGE
DEP_INCLUDES  += -I${ADAGE_HOME}/include
DEP_LIBRARIES += -L${ADAGE_HOME}/lib -lAdage
DEP_TARGETS   += ${ADAGE_HOME}/lib/libAdage.so


# epic-analysis objects
#################################################################

PREFIX                  := ${EPIC_ANALYSIS_HOME}/lib
EPIC_ANALYSIS_LIB       := $(PREFIX)/lib$(PROJECT).so
EPIC_ANALYSIS_PCM       := $(PROJECT)Dict_rdict.pcm
EPIC_ANALYSIS_DICT      := src/$(PROJECT)Dict.cxx
EPIC_ANALYSIS_LINKDEF   := src/LinkDef.h
EPIC_ANALYSIS_INCLUDES  := -I${EPIC_ANALYSIS_HOME}/src
EPIC_ANALYSIS_LIBRARIES := -L$(PREFIX) -l$(PROJECT)

# epic-analysis source code (with $(EPIC_ANALYSIS_DICT) and $(EPIC_ANALYSIS_LINKDEF) moved to end of lists for rootcling)
#################################################################

# move $(EPIC_ANALYSIS_DICT) and $(EPIC_ANALYSIS_LINKDEF) to end of lists for rootcling
SOURCES := $(filter-out $(EPIC_ANALYSIS_DICT),    $(wildcard src/*.cxx) $(wildcard src/sfset/*.cxx) $(wildcard src/interp/*.cxx)) $(EPIC_ANALYSIS_DICT)
HEADERS := $(filter-out $(EPIC_ANALYSIS_LINKDEF), $(wildcard src/*.h)   $(wildcard src/sfset/*.h)   $(wildcard src/interp/*.h)    $(wildcard src/interp/*.ipp)) $(EPIC_ANALYSIS_LINKDEF)

# filter out Delphes-dependent code (if needed)
ifdef EXCLUDE_DELPHES
	SOURCES := $(filter-out src/AnalysisDelphes.cxx, $(SOURCES))
	HEADERS := $(filter-out src/AnalysisDelphes.h,   $(HEADERS))
	ROOTCLING += -D=EXCLUDE_DELPHES
endif


# epic-analysis builds
#################################################################

all: epic-analysis
all-clean: deps-clean clean

debug: FLAGS += -g
debug: DEP_RECIPE = debug
debug: clean all
release: FLAGS += -O3
release: DEP_RECIPE = release
release: clean all

epic-analysis: deps epic-analysis-header $(EPIC_ANALYSIS_LIB) hpc
	@echo "Done.\n"
epic-analysis-header:
	@echo "\n===== $(PROJECT) ====="

$(EPIC_ANALYSIS_LIB): $(EPIC_ANALYSIS_DICT) $(HEADERS) $(SOURCES)
	@echo "\n----- $(PROJECT) library -----"
	@if [ -z "${EXCLUDE_DELPHES}" ]; then \
		ln -svf ${DELPHES_HOME}/external ./; \
	fi
	@echo "target: $@"
	$(CXX) $(SOURCES) -shared -o $@ $(FLAGS) $(DEP_LIBRARIES) $(DEP_INCLUDES) -I${EPIC_ANALYSIS_HOME}

$(EPIC_ANALYSIS_DICT): $(HEADERS) $(DEP_TARGETS)
	@echo "\n----- $(PROJECT) dictionary -----"
	@echo "target: $@"
	$(ROOTCLING) -f $@ $(DEP_INCLUDES) $(HEADERS)
	@mkdir -p $(PREFIX)
	@mv src/$(EPIC_ANALYSIS_PCM) $(PREFIX)/

# dependency builds
#################################################################

deps: delphes mstwpdf adage
deps-clean: delphes-clean mstwpdf-clean adage-clean
delphes:
	@if [ -z "${EXCLUDE_DELPHES}" ]; then \
		echo "\n===== $@ ====="; \
		$(MAKE) -C ${DELPHES_HOME}; \
	fi
delphes-clean:
	@if [ -z "${EXCLUDE_DELPHES}" ]; then \
		echo "\n===== $@ ====="; \
		$(MAKE) -C ${DELPHES_HOME} clean; \
	fi
mstwpdf:
	@echo "\n===== $@ ====="
	$(MAKE) -C ${MSTWPDF_HOME} $(DEP_RECIPE)
mstwpdf-clean:
	@echo "\n===== $@ ====="
	$(MAKE) -C ${MSTWPDF_HOME} clean
adage:
	@echo "\n===== $@ ====="
	$(MAKE) -C ${ADAGE_HOME} $(DEP_RECIPE)
adage-clean:
	@echo "\n===== $@ ====="
	$(MAKE) -C ${ADAGE_HOME} clean


# HPC executables
#################################################################

HPC_SOURCES     := $(basename $(wildcard hpc/src/*.cpp))
HPC_EXECUTABLES := $(addsuffix .exe, $(HPC_SOURCES))

hpc: hpc-header $(HPC_EXECUTABLES)
hpc-header:
	@echo "\n===== HPC executables ====="
hpc/src/%.exe: hpc/src/%.cpp $(EPIC_ANALYSIS_LIB) 
	@echo "----- build $@.o -----"
	$(CXX) -c $< -o $@.o $(FLAGS) $(DEP_INCLUDES) $(EPIC_ANALYSIS_INCLUDES)
	@echo "--- make executable $@"
	$(CXX) -o $@ $@.o $(DEP_LIBRARIES) $(EPIC_ANALYSIS_LIBRARIES)
	$(RM) $@.o


# clean
#################################################################
clean:
	@echo "\n===== CLEAN $(PROJECT) ====="
	$(RM) $(EPIC_ANALYSIS_LIB) $(PREFIX)/$(EPIC_ANALYSIS_PCM) $(EPIC_ANALYSIS_DICT) $(HPC_EXECUTABLES)
