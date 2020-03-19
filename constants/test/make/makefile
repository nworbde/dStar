include $(MESA_DIR)/utils/makefile_header

DSTAR_LIB_DIR = ../../../lib

LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lconstants
LOAD_OTHER = -L$(DSTAR_LIB_DIR) -lexceptions -L$(MESA_LIB_DIR) -lconst -lutils

TEST = tester
TEST_OBJS = test_constants.o

default: install

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) $<

install: $(TEST)
	cp $(TEST) ..

clean:
	-@rm -f *.o *.mod $(TEST)
