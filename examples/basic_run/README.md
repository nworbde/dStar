INSTRUCTIONS FOR USING THIS EXAMPLE
===================================

This directory contains an example run of an accretion outburst/quiescent cooling. 

1.  This directory may be copied anywhere on your computer. Note, however, that you need write permission to dStar/data to run as the code caches datafiles there.  I'll change this in the near future to allow you to set the cache location (useful for running on multi-user machines).

2.  Ensure that the environment variable `DSTAR_DIR` is set. This is used in the makefile,

        DSTAR_LIB_DIR = $(DSTAR_DIR)/lib
        DSTAR_INC = $(DSTAR_DIR)/include

3.  Change the parameters as desired in the `inlist`.
    
4.  Build the code: `./mk`
    
5.  And run it: `./run_dStar -D<path to dStar> -I<inlist name>`
    The command now takes two optional arguments, the path to the root directory of dStar and the name of the inlist file.  You no longer need to edit the src code for that.  The path to dStar is normally read from the environment variable `DSTAR_DIR`, but can be overwritten here. The default value for the inlist is
    
        default_inlist_file = 'inlist'
    
Happy modeling!
