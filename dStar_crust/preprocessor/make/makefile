include $(MESA_DIR)/utils/makefile_header

DSTAR_LIB_DIR = ../../../lib
LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -ldstar_crust
LOAD_OTHER = -L$(DSTAR_LIB_DIR) -ldstar_eos -lsuperfluid -lnucchem -lconstants -lexceptions -L$(MESA_LIB_DIR) $(LOAD_MESA_NUMERICS)
DSTAR_INC = ../../../include

HZ90_GENERATOR = generate_HZ90
HZ90_GENERATOR_OBJS =  HZ90_comp.o composition_handler.o generate_HZ90.o

PREPROCESSOR = process_abuntime
PREPROCESSOR_OBJS = HZ90_comp.o abuntime.o composition_handler.o process_abuntime.o

PRINTER = print_composition
PRINTER_OBJS = composition_handler.o print_composition.o

default: install

$(PREPROCESSOR) : $(PREPROCESSOR_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(PREPROCESSOR_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

$(PRINTER) : $(PRINTER_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(PRINTER_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

$(HZ90_GENERATOR): $(HZ90_GENERATOR_OBJS)
	$(LOADER) $(FCopenmp) -o $@ $(HZ90_GENERATOR_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) -I$(DSTAR_INC) $<

install: $(PREPROCESSOR) $(PRINTER) $(HZ90_GENERATOR)
	cp $(PREPROCESSOR) $(PRINTER) $(HZ90_GENERATOR) ..

clean:
	-@rm -f *.o *.mod $(PREPROCESSOR) $(PRINTER) $(HZ90_GENERATOR)
