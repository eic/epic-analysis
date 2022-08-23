# compiler and flags
CXX = g++
FLAGS = -g -Wno-deprecated -fPIC -fno-inline -Wno-write-strings
FLAGS += -fmax-errors=3
# FLAGS += -O0 # TODO: use optimization?

# ROOT
ROOT_INCLUDES := -I$(shell root-config --incdir)
ROOT_LIBS     := $(shell root-config --glibs)

#-----------------------------------------------

# dirs
SRC_DIR := src
INC_DIR := include/adage
LIB_DIR := lib

# targets
ADAGE   := Adage
LIB     := $(LIB_DIR)/lib$(ADAGE).so
DICT    := $(SRC_DIR)/$(ADAGE)Dict.cxx
PCM     := $(ADAGE)Dict_rdict.pcm
LINKDEF := include/LinkDef.h

# sources
HEADERS := $(wildcard $(INC_DIR)/*.h)
SOURCES := $(filter-out $(DICT), $(wildcard $(SRC_DIR)/*.cxx)) $(DICT)

#-----------------------------------------------

$(LIB): $(SOURCES) $(HEADERS) $(LINKDEF)
	@echo ">>>>> build $(ADAGE) library"
	@echo "      $(CURDIR)/$@"
	$(CXX) $(SOURCES) -shared -o $@ $(FLAGS) -I$(INC_DIR) $(ROOT_INCLUDES) $(ROOT_LIBS)

$(DICT): $(HEADERS) $(LINKDEF)
	@echo ">>>>> generate $(ADAGE) dictionary"
	@echo "      $(CURDIR)/$@"
	@mkdir -p $(LIB_DIR)
	@cd $(INC_DIR) && rootcling -f $(CURDIR)/$@ $(ROOT_INCLUDES) $(HEADERS:$(INC_DIR)/%=%) $(CURDIR)/$(LINKDEF)
	@mv $(SRC_DIR)/$(PCM) $(LIB_DIR)

clean:
	@echo ">>>>> clean $(ADAGE)"
	$(RM) $(LIB) $(DICT) $(LIB_DIR)/$(PCM)
