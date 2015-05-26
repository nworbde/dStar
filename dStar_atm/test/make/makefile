include $(MESA_DIR)/utils/makefile_header

DSTAR_INC=../../../include
DSTAR_LIB_DIR=../../../lib

LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -ldstar_atm

LOAD_OTHER = -L$(DSTAR_LIB_DIR) -lconductivity -lneutrino -lconstants -ldstar_eos -lsuperfluid  -lnucchem -L$(MESA_LIB_DIR) $(LOAD_MESA_NUMERICS)

TEST = tester
TEST_OBJS = test_atm.o #test_bc09.o

default: install

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) -I$(DSTAR_INC) $<

install: $(TEST)
	cp $(TEST) ..

clean:
	-@rm -f *.o *.mod $(TEST)
