include config.mk

# sidis-eic targets
all: deps sidis-eic
all-clean: deps-clean clean
sidis-eic:
	ln -sf ${DELPHES_HOME}/external ./
	@cd src; make
clean:
	@cd src; make clean

# dependency targets
deps: delphes mstwpdf
deps-clean: delphes-clean mstwpdf-clean
delphes:
	@cd ${DELPHES_HOME}; make
delphes-clean:
	@cd ${DELPHES_HOME}; make clean
mstwpdf:
	@cd ${MSTWPDF_HOME}; make
mstwpdf-clean:
	@cd ${MSTWPDF_HOME}; make clean
