INSTRUCTIONS FOR USING THIS EXAMPLE
===================================

This directory contains a demonstration of how the thermal conductivity and heating profile affects the quiescent cooling. This material was used in a talk I gave during the INT workshop 16-2b, "The Phases of Dense Matter", 20 July 2016.

1.  This directory may be copied anywhere on your computer. Note, however, that you need write permission to dStar/data to run as the code caches datafiles there.  I'll change this in the near future to allow you to set the cache location (useful for running on multi-user machines).  If you do copy this directory to another location on the computer, you will need to make the following edits.

    1.  Edit make/makefile.  Replace the following macros with the correct paths

            DSTAR_LIB_DIR = ../../dStar/lib
            DSTAR_INC = ../../dStar/include

    2.  Edit the line

            DSTAR_DIR='../../'

        in both `mk` and `rn` to point to the dStar directory.
    
2.  Build the code: `./mk`
    
3.  And run it: `./rn`

4.  You can view the outputs by opening the Jupyter notebook `INT2016-crusts.ipynb`

