include config.mk

# dependencies
DEPS += -Isrc
#LIBS += -L. -l$(SIDIS_EIC)

# sidis-eic targets
sidis-eic:
	ln -sf ${DELPHES_HOME}/external ./
	@cd ${MSTWPDF_HOME}; make
	@cd src; make
clean:
	@cd ${MSTWPDF_HOME}; make clean
	@cd src; make clean
all:
	make deps
	make sidis-eic
all-clean:
	make deps-clean
	make clean

# dependency targets
deps:
	make delphes
deps-clean:
	make delphes-clean
delphes:
	@cd ${DELPHES_HOME}; make
delphes-clean:
	@cd ${DELPHES_HOME}; make clean
