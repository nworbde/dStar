include $(MESA_DIR)/utils/makefile_header

DSTAR_LIB = ../../../lib
DSTAR_INC = ../../../include

LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lchem 
LOAD_OTHER = -L$(DSTAR_LIB) -lnetJina -L$(MESA_LIB_DIR) -lconst -lutils

TEST = tester
TEST_OBJS = test_chem.o

default: install

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) -I$(DSTAR_INC) $<

install: $(TEST)
	cp $(TEST) ..

clean:
	-@rm -f *.o *.mod $(TEST)
