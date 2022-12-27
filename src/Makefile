include ../config.mk

EPIC_ANALYSIS = EpicAnalysis
ROOTCLING = rootcling

# targets
PREFIX  := ${EPIC_ANALYSIS_HOME}/lib
LIB     := lib$(EPIC_ANALYSIS).so
PCM     := $(EPIC_ANALYSIS)Dict_rdict.pcm
DICT    := $(EPIC_ANALYSIS)Dict.cxx
LINKDEF := LinkDef.h

# source code (with $(DICT) and $(LINKDEF) moved to end of lists for rootcling)
SOURCES := $(filter-out $(DICT),    $(wildcard *.cxx) $(wildcard sfset/*.cxx) $(wildcard interp/*.cxx)) $(DICT)
HEADERS := $(filter-out $(LINKDEF), $(wildcard *.h)   $(wildcard sfset/*.h)   $(wildcard interp/*.h)    $(wildcard interp/*.ipp)) $(LINKDEF)
ifdef EXCLUDE_DELPHES
	SOURCES := $(filter-out AnalysisDelphes.cxx, $(SOURCES))
	HEADERS := $(filter-out AnalysisDelphes.h,   $(HEADERS))
	ROOTCLING += -D=EXCLUDE_DELPHES
endif

#-----------------------------------------------

$(PREFIX)/$(LIB): $(DICT) $(HEADERS) $(SOURCES)
	@echo ">>>>> build $(EPIC_ANALYSIS) library"
	@echo "      $@"
	$(CXX) $(SOURCES) -shared -o $@ $(FLAGS) $(DEPS) $(LIBS)

$(DICT): $(HEADERS)
	@echo ">>>>> generate $(EPIC_ANALYSIS) dictionary"
	@echo "      $@"
	$(ROOTCLING) -f $@ $(DEPS) $^
	@mkdir -p $(PREFIX)
	@mv $(PCM) $(PREFIX)/

clean:
	$(RM) $(PREFIX)/$(LIB) $(PREFIX)/$(PCM) $(DICT)
