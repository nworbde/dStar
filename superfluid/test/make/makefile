include $(MESA_DIR)/utils/makefile_header

LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lsuperfluid
LOAD_OTHER = -L$(MESA_LIB_DIR) -linterp_1d -lconst -lutils

TEST = tester
TEST_OBJS = test_sf.o

default: install

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) $<

install: $(TEST)
	cp $(TEST) ..

clean:
	-@rm -f *.o *.mod $(TEST)
