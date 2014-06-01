include $(MESA_DIR)/utils/makefile_header

LIB = libcrust_eos.a
LIB_OBJS = dStar_eos_def.o funct_fermi1.o neutron.o helm_alloc.o helm_core.o electron.o ion.o crust_eos_lib.o

VPATH = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR)
DSTAR_INC=../../include
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
	cp $(LIB) $(BASE_DIR)/lib ; cp *.mod $(BASE_DIR)/include
	cd $(BASE_DIR)/lib ; ranlib $(LIB)

clean:
	-@rm -f *.o *.mod $(LIB)
