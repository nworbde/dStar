include $(MESA_DIR)/utils/makefile_header

DSTAR_LIB_DIR = /path/to/dStar/lib
LOAD_OTHER = -L$(DSTAR_LIB_DIR) -ldstar_crust -lconductivity -lneutrino -ldstar_eos -lsuperfluid -ldstar_atm -lnucchem -lconstants -lexceptions -L$(MESA_LIB_DIR) $(LOAD_MESA_NUMERIC)
DSTAR_INC = /path/to/dStar/include

EXEC = run_MRcurve
EXEC_OBJS = run.o

default: install

$(EXEC) : $(EXEC_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(EXEC_OBJS) $(LOAD_OTHER)

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) -I$(DSTAR_INC) $<

install: $(EXEC)
	cp $(EXEC) ..
	mkdir -p ../LOGS

clean:
	-@rm -f *.o *.mod $(EXEC)
