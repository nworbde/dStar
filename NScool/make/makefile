include $(MESA_DIR)/utils/makefile_header

LIB = libnscool.a
LIB_OBJS = NScool_def.o NScool_private_def.o NScool_storage.o NScool_crust_tov.o NScool_ctrls_io.o NScool_create_model.o  NScool_init.o NScool_profile.o NScool_history.o NScool_terminal.o NScool_evolve.o NScool_lib.o argparse.o

MOD_JOB_DIR=../job
VPATH = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR):$(MOD_JOB_DIR)
DSTAR_INC=../../include
DSTAR_LIB=../../lib
COMPILE = $(COMPILE_TO_TEST) $(FCfree) -I$(DSTAR_INC) -I../defaults
COMPILE_FIXED = $(COMPILE_TO_TEST) $(FCfixed) -I$(DSTAR_INC)

%.o:%.f
	$(COMPILE) $<

$(LIB): $(LIB_OBJS)
	$(LIB_TOOL) $(LIB) $(LIB_OBJS)

install: $(LIB)
	cp $(LIB) $(DSTAR_LIB) ; cp *.mod $(DSTAR_INC)
	cd $(DSTAR_LIB) ; ranlib $(LIB)

clean:
	-@rm -f *.o *.mod $(LIB)
