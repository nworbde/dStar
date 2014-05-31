include $(MESA_DIR)/utils/makefile_header

DSTAR_LIB_DIR = ../../../lib
LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lneutrino 
LOAD_OTHER = -L$(DSTAR_LIB_DIR) -lsuperfluid -lnucchem -lnetJina -lconstants -L$(MESA_LIB_DIR) -linterp_1d -lconst -lutils
DSTAR_INC = ../../../include

TEST = tester
TEST_OBJS = test_neutrino.o

default: install

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

test_crust : test_crust.o
	$(LOADER) $(FCopenmp) -o $@ test_crust.o $(LOAD_LOCAL) $(LOAD_OTHER)

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) -I$(DSTAR_INC) $<

install: $(TEST) test_crust
	cp $(TEST) test_crust ..

clean:
	-@rm -f *.o *.mod $(TEST)
