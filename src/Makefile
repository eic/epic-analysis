include ../config.mk

# rootcling dictionary
DICT := $(SIDIS-EIC)Dict.cxx
PCM := $(SIDIS-EIC)Dict_rdict.pcm
LINKDEF = LinkDef.h


# source code (with $(DICT) and $(LINKDEF) moved to end of lists for rootcling)
SOURCES := $(filter-out $(DICT), $(wildcard *.cxx) $(wildcard sfset/*.cxx) $(wildcard interp/*.cxx)) $(DICT)
HEADERS := $(filter-out $(LINKDEF), $(wildcard *.h) $(wildcard sfset/*.h) $(wildcard interp/*.h) $(wildcard interp/*.ipp)) $(LINKDEF)


#-----------------------------------------------


$(SIDIS-EIC-OBJ): $(DICT) $(HEADERS) $(SOURCES)
	@echo "----- build $@ -----"
	$(CXX) $(SOURCES) -shared -o ../$@ $(FLAGS) $(DEPS) $(LIBS)

$(DICT): $(HEADERS)
	@echo "----- generate dictionary $@ -----"
	@rootcling -f $@ $(DEPS) $^
	mv $(PCM) ../

clean:
	$(RM) ../$(SIDIS-EIC-OBJ) $(DICT) ../$(PCM)
