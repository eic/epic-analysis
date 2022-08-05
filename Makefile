include config.mk

# dependencies
DEPS += -Isrc
#LIBS += -L. -l$(SIDIS-EIC)

# sidis-eic targets
sidis-eic:
	@cd mstwpdf; make
	@cd src; make
clean:
	@cd mstwpdf; make clean
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
