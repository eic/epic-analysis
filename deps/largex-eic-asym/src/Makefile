include ../config.mk

# rootcling dictionary
DICT := $(LARGEXASYM)Dict.cxx
PCM := $(LARGEXASYM)Dict_rdict.pcm
LINKDEF = LinkDef.h


# source code (with $(DICT) and $(LINKDEF) moved to end of lists for rootcling)
SOURCES := $(filter-out $(DICT), $(wildcard *.cxx)) $(DICT)
HEADERS := $(filter-out $(LINKDEF), $(wildcard *.h)) $(LINKDEF)


#-----------------------------------------------


$(LARGEXASYMOBJ): $(DICT) $(HEADERS) $(SOURCES)
	@echo "----- build $@ -----"
	$(CXX) $(SOURCES) -shared -o ../$@ $(FLAGS) $(DEPS) $(LIBS)

$(DICT): $(HEADERS)
	@echo "----- generate dictionary $@ -----"
	@echo "WARNING: -Wno-defaulted-function-deleted is used to silence brufit warnings"
	@rootcling -f $@ -Wno-defaulted-function-deleted $(DEPS) $^
	mv $(PCM) ../

clean:
	$(RM) ../$(LARGEXASYMOBJ) $(DICT) ../$(PCM)
