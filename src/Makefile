include ../config.mk

SIDIS_EIC = Sidis-eic
ROOTCLING = rootcling

# targets
PREFIX  := ${SIDIS_EIC_HOME}/lib
LIB     := lib$(SIDIS_EIC).so
PCM     := $(SIDIS_EIC)Dict_rdict.pcm
DICT    := $(SIDIS_EIC)Dict.cxx
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
	@echo ">>>>> build $(SIDIS_EIC) library"
	@echo "      $@"
	$(CXX) $(SOURCES) -shared -o $@ $(FLAGS) $(DEPS) $(LIBS)

$(DICT): $(HEADERS)
	@echo ">>>>> generate $(SIDIS_EIC) dictionary"
	@echo "      $@"
	$(ROOTCLING) -f $@ $(DEPS) $^
	@mkdir -p $(PREFIX)
	@mv $(PCM) $(PREFIX)/

clean:
	$(RM) $(PREFIX)/$(LIB) $(PREFIX)/$(PCM) $(DICT)
