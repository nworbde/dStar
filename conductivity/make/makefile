include $(MESA_DIR)/utils/makefile_header

LIB = libconductivity.a
LIB_OBJS = conductivity_def.o eval_conductivity.o eval_core_conductivity.o conductivity_lib.o
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
