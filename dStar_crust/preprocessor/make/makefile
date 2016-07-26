include $(MESA_DIR)/utils/makefile_header

DSTAR_LIB_DIR = ../../../lib
LOAD_OTHER = -L$(DSTAR_LIB_DIR) -ldstar_eos -lsuperfluid -lnucchem -lconstants -L$(MESA_LIB_DIR) $(LOAD_MESA_NUMERICS)
DSTAR_INC = ../../../include

PREPROCESSOR = process_abuntime
PREPROCESSOR_OBJS = abuntime.o process_abuntime.o

default: install

$(PREPROCESSOR) : $(PREPROCESSOR_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(PREPROCESSOR_OBJS) $(LOAD_OTHER)

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) -I$(DSTAR_INC) $<

install: $(PREPROCESSOR)
	cp $(PREPROCESSOR) ..

clean:
	-@rm -f *.o *.mod $(PREPROCESSOR)
