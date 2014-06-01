include $(MESA_DIR)/utils/makefile_header

LIB = libdstar_eos.a
LIB_OBJS = 	dStar_eos_def.o dStar_eos_private_def.o funct_fermi1.o \
			neutron.o helm_alloc.o helm_core.o electron.o ion.o dStar_eos_lib.o

VPATH = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR)
DSTAR_INC=../../include
DSTAR_LIB=../../lib
COMPILE = $(COMPILE_TO_TEST) $(FCfree) -I$(DSTAR_INC)
COMPILE_FIXED = $(COMPILE_TO_TEST) $(FCfixed) -I$(DSTAR_INC)

%.o:%.f
	$(COMPILE) $<

$(LIB): $(LIB_OBJS)
	$(LIB_TOOL) $(LIB) $(LIB_OBJS)

helm_alloc.o: $(MOD_PRIVATE_DIR)/helm_alloc.f
	 $(COMPILE_FIXED) $<

helm_core.o: $(MOD_PRIVATE_DIR)/helm_core.f
	$(COMPILE_FIXED) $<

generate_acprs_table: acprs.o
	$(FC) -o $@ acprs.o -L. -lcrust_eos -L../../lib -L../../mesa_numerics/lib -lchem -lalert -lsuperfluid -linterp_1d

install: $(LIB)
	cp $(LIB) $(DSTAR_LIB) ; cp *.mod $(DSTAR_INC)
	cd $(DSTAR_LIB) ; ranlib $(LIB)

clean:
	-@rm -f *.o *.mod $(LIB)
