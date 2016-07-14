include $(MESA_DIR)/utils/makefile_header

LIB = libdstar_crust.a
LIB_OBJS = 	dStar_crust_def.o HZ90.o abuntime.o dStar_crust_mod.o dStar_crust_lib.o
DSTAR_INC=../../include
DSTAR_LIB=../../lib

VPATH = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR)

COMPILE = $(COMPILE_TO_TEST) $(FCfree) -I$(DSTAR_INC)

%.o:%.f
	$(COMPILE) $<

default: $(LIB)

$(LIB): $(LIB_OBJS)
	$(LIB_TOOL) $(LIB) $(LIB_OBJS)

install: $(LIB)
	cp $(LIB) $(DSTAR_LIB) ; cp *.mod $(DSTAR_INC)
	cd $(DSTAR_LIB) ; ranlib $(LIB)

clean:
	-@rm -f *.o *.mod $(LIB)
