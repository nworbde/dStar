include $(MESA_DIR)/utils/makefile_header

DSTAR_LIB_DIR = ../../../lib
DSTAR_INC = ../../../include
LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lsuperfluid
LOAD_OTHER = -L$(DSTAR_LIB_DIR) -lconstants -lexceptions -L$(MESA_LIB_DIR) -linterp_1d -lconst -lutils

TEST = tester
TEST_OBJS = test_sf.o

default: install

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) -I$(DSTAR_INC) $<

install: $(TEST)
	cp $(TEST) ..

clean:
	-@rm -f *.o *.mod $(TEST)
