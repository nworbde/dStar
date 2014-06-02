include $(MESA_DIR)/utils/makefile_header

DSTAR_INC=../../../include
DSTAR_LIB=../../../lib

LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -ldstar_atm
LOAD_OTHER = -L$(DSTAR_LIB) -lconstants -L$(MESA_LIB_DIR) -linterp_1d -lconst -lutils

TEST = tester
TEST_OBJS = test_atm.o

default: install

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) -I$(DSTAR_INC) $<

install: $(TEST)
	cp $(TEST) ..

clean:
	-@rm -f *.o *.mod $(TEST)
