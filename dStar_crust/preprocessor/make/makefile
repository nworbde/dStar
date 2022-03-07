include $(MESA_DIR)/utils/makefile_header

DSTAR_LIB_DIR = ../../../lib
LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -ldstar_crust
LOAD_OTHER = -L$(DSTAR_LIB_DIR) -ldstar_eos -lsuperfluid -lnucchem -lconstants -lexceptions -L$(MESA_LIB_DIR) $(LOAD_MESA_NUMERICS)
DSTAR_INC = ../../../include

HZ90_GENERATOR = generate_HZ90
HZ90_GENERATOR_OBJS =  HZ90_comp.o generate_HZ90.o

PREPROCESSOR = process_abuntime
PREPROCESSOR_OBJS = HZ90_comp.o abuntime.o process_abuntime.o

PRINTER = print_composition
PRINTER_OBJS = print_composition.o

ANALYZER = analyze_crust
ANALYZER_OBJS = analyze_crust.o

default: install

$(HZ90_GENERATOR): $(HZ90_GENERATOR_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(HZ90_GENERATOR_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

$(PREPROCESSOR) : $(PREPROCESSOR_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(PREPROCESSOR_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

$(PRINTER) : $(PRINTER_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(PRINTER_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

$(ANALYZER) : $(ANALYZER_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(ANALYZER_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) -I$(DSTAR_INC) $<

install: $(HZ90_GENERATOR) $(PREPROCESSOR) $(PRINTER) $(ANALYZER)
	cp $(HZ90_GENERATOR) $(PREPROCESSOR) $(PRINTER) $(ANALYZER) ..

clean:
	-@rm -f *.o *.mod $(HZ90_GENERATOR) $(PREPROCESSOR) $(PRINTER) $(ANALYZER)
