include $(MESA_DIR)/utils/makefile_header

DSTAR_LIB_DIR = ../../../lib
LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -ldstar_eos
LOAD_OTHER = -L$(DSTAR_LIB_DIR) -lsuperfluid -lnucchem  -lconstants -L$(MESA_LIB_DIR) -linterp_1d -lconst -lutils
DSTAR_INC = ../../../include

TEST = tester
TEST_OBJS = test_crust_eos.o

CORE_TEST = core_test
CORE_TEST_OBJS = test_core_eos.o

default: install

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

$(CORE_TEST) : $(CORE_TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(CORE_TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) -I$(DSTAR_INC) $<

install: $(TEST) $(CORE_TEST)
	cp $(TEST) $(CORE_TEST) ..

clean:
	-@rm -f *.o *.mod $(TEST) $(CORE_TEST)
