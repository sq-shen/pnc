all: bicm_pnc

bicm_pnc:
	cd bicm_pnc && $(MAKE)

clean:
	cd bicm_pnc && $(MAKE) clean


.PHONY: clean bicm_pnc
