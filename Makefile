PKG_DIR = $(PWD)
SUB_DIR = $(PKG_DIR)/submodule
AFT_DIR = $(SUB_DIR)/AfterPhylo
BFR_DIR = $(SUB_DIR)/BeforePhylo

all: modules

modules:
	git submodule update --init --recursive
	cp $(AFT_DIR)/AfterPhylo.pl ./bin/; chmod 755 ./bin/AfterPhylo.pl
	cp $(BFR_DIR)/BeforePhylo.pl ./bin/; chmod 755 ./bin/BeforePhylo.pl

clean:
	-rm -f ./bin/AfterPhylo.pl
	-rm -f ./bin/BeforePhylo.pl
