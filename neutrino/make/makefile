include $(MESA_DIR)/utils/makefile_header

LIB = libneutrino.a
LIB_OBJS = 	neutrino_def.o \
			eval_core_neutrino.o \
			eval_crust_neutrino.o \
			neutrino_lib.o

DSTAR_LIB = ../../lib
DSTAR_INC = ../../include

VPATH = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR)

$(LIB): $(LIB_OBJS)
	$(LIB_TOOL) $(LIB) $(LIB_OBJS)

COMPILE = $(COMPILE_TO_TEST) $(FCfree) -I$(DSTAR_INC)

%.o:%.f
	$(COMPILE) $<

default: $(LIB)

install: $(LIB)
	cp $(LIB) $(DSTAR_LIB) ; cp *.mod $(DSTAR_INC)
	cd $(DSTAR_LIB) ; ranlib $(LIB)

clean:
	-@rm -f *.o *.mod $(LIB)
