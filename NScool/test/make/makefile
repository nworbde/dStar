include $(MESA_DIR)/utils/makefile_header

DSTAR_LIB_DIR = ../../../lib
LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lnscool
LOAD_OTHER = -L$(DSTAR_LIB_DIR) -ldstar_crust -lsuperfluid -lconductivity -lneutrino -ldstar_eos -ldstar_atm -lnucchem -lconstants -L$(MESA_LIB_DIR) $(LOAD_MESA_NUMERICS)
DSTAR_INC = ../../../include

TEST = tester
TEST_OBJS = test_NScool.o #

default: install

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) -I$(DSTAR_INC) $<

install: $(TEST)
	cp $(TEST) ..
	mkdir -p ../LOGS

clean:
	-@rm -f *.o *.mod $(TEST)
