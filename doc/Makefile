LANGUAGE=DE

all: install
install: install_theory install_doxygen

install_theory:
	$(MAKE) -C theory/ install LANG=$(LANGUAGE)

install_doxygen:
	$(MAKE) -C doxygen/ install
