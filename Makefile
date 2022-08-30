include config.mk

# sidis-eic targets
all: deps sidis-eic
all-clean: deps-clean clean
sidis-eic:
	@echo "\n===== SIDIS-EIC ====="
	@if [ -z "${EXCLUDE_DELPHES}" ]; then \
		ln -svf ${DELPHES_HOME}/external ./; \
	fi
	@cd src; $(MAKE)
clean:
	@echo "\n===== CLEAN SIDIS-EIC ====="
	@cd src; $(MAKE) clean

# dependency targets
deps: delphes mstwpdf adage
deps-clean: delphes-clean mstwpdf-clean adage-clean
delphes:
	@if [ -z "${EXCLUDE_DELPHES}" ]; then \
		echo "\n===== DELPHES ====="; \
		cd ${DELPHES_HOME}; $(MAKE); \
	fi
delphes-clean:
	@if [ -z "${EXCLUDE_DELPHES}" ]; then \
		echo "\n===== CLEAN DELPHES ====="; \
		cd ${DELPHES_HOME}; $(MAKE) clean; \
	fi
mstwpdf:
	@echo "\n===== MSTWPDF ====="
	@cd ${MSTWPDF_HOME}; $(MAKE)
mstwpdf-clean:
	@echo "\n===== CLEAN MSTWPDF ====="
	@cd ${MSTWPDF_HOME}; $(MAKE) clean
adage:
	@echo "\n===== ADAGE ====="
	@cd ${ADAGE_HOME}; $(MAKE)
adage-clean:
	@echo "\n===== CLEAN ADAGE ====="
	@cd ${ADAGE_HOME}; $(MAKE) clean
