include $(MESA_DIR)/utils/makefile_header

LIB = libmrcurve.a
LIB_OBJS = dStar_core_def.o brown_skyrme.o load_skyrme_parameters.o \
			dStar_core_mod.o dStar_core_lib.o dStar_core_tov.o

VPATH = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR)
DSTAR_INC=../../include
DSTAR_LIB=../../lib
COMPILE = $(COMPILE_TO_TEST) $(FCfree) -I$(DSTAR_INC)

%.o:%.f
	$(COMPILE) $<

$(LIB): $(LIB_OBJS)
	$(LIB_TOOL) $(LIB) $(LIB_OBJS)

install: $(LIB)
	cp $(LIB) $(DSTAR_LIB) ; cp *.mod $(DSTAR_INC)
	cd $(DSTAR_LIB) ; ranlib $(LIB)

clean:
	-@rm -f *.o *.mod $(LIB)
