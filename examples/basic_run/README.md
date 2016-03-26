INSTRUCTIONS FOR USING THIS EXAMPLE
===================================

This directory contains an example run of an accretion outburst/quiescent cooling. 

1. This directory may be copied anywhere on your computer. Note, however, that you need write permission to dStar/data to run as the code caches datafiles there.  I'll change this in the near future to allow you to set the cache location (useful for running on multi-user machines).

2. Edit make/makefile.  Replace the following macros with the correct paths

        DSTAR_LIB_DIR = /path/to/dStar/lib
        DSTAR_INC = /path/to/dStar/include

3. Edit src/run.f.  Replace this variable with the correct path

        character(len=*), parameter :: my_dStar_dir = '/path/to/local/dStar'

4. Change the parameters as desired in the `inlist`.
    
5. Build the code: `./mk`
    
6. And run it: `./run_dStar`
    
Happy modeling!
