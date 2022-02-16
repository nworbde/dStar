include $(MESA_DIR)/utils/makefile_header

DSTAR_LIB_DIR = ../../../lib
LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lconductivity
LOAD_OTHER = -L$(DSTAR_LIB_DIR) -lneutrino -lconstants -ldstar_eos -lsuperfluid  -lnucchem -lexceptions -L$(MESA_LIB_DIR) -linterp_1d -lconst -lutils $(LOAD_MESA_NUMERICS)
DSTAR_INC = ../../../include

TEST = tester
TEST_OBJS = test_cond.o

EXPERIMENTAL = tabulate_conductivity test_tabulated_mixture tabulate_mixture

default: install

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)
	
tabulate_conductivity: tabulate_conductivity.o
	$(LOADER) $(FCopenmp) -o $@ $< $(LOAD_LOCAL) $(LOAD_OTHER)

test_tabulated_mixture: test_tabulated_mixture.o
	$(LOADER) $(FCopenmp) -o $@ $< $(LOAD_LOCAL) $(LOAD_OTHER)

tabulate_mixture: tabulate_mixture.o
	$(LOADER) $(FCopenmp) -o $@ $< $(LOAD_LOCAL) $(LOAD_OTHER)

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) -I$(DSTAR_INC) $<

install: $(TEST) $(EXPERIMENTAL)
	cp $(TEST) $(EXPERIMENTAL) ..

clean:
	-@rm -f *.o *.mod $(TEST) $(EXPERIMENTAL)
