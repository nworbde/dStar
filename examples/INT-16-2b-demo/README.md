INSTRUCTIONS FOR USING THIS EXAMPLE
===================================

This directory contains a demonstration of how the thermal conductivity and heating profile affects the quiescent cooling. This material was used in a talk I gave during the INT workshop 16-2b, "The Phases of Dense Matter", 20 July 2016.

1.  This directory may be copied anywhere on your computer. Note, however, that you need write permission to dStar/data to run as the code caches datafiles there.  I'll change this in the near future to allow you to set the cache location (useful for running on multi-user machines).  If you do copy this directory to another location on the computer, you will need to make the following edits.

    Ensure that the environment variable `DSTAR_DIR` is set. This is used in the makefile,

        DSTAR_LIB_DIR = $(DSTAR_DIR)/lib
        DSTAR_INC = $(DSTAR_DIR)/include

    
2.  Build the code: `./mk`
    
3.  And run it: `./rn`

4.  You can view the outputs by opening the Jupyter notebook `INT2016-crusts.ipynb`

