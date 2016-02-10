include $(MESA_DIR)/utils/makefile_header

DSTAR_LIB_DIR = ../../../lib
LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -ldstar_core
LOAD_OTHER = -L$(DSTAR_LIB_DIR) -ldstar_eos -lsuperfluid -lnucchem -lconstants -L$(MESA_LIB_DIR) $(LOAD_MESA_NUMERICS)
DSTAR_INC = ../../../include

TEST = tester
TEST_OBJS = test_core.o

default: install

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) -I$(DSTAR_INC) $<

install: $(TEST)
	cp $(TEST) ..

clean:
	-@rm -f *.o *.mod $(TEST)
