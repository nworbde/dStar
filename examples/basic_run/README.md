INSTRUCTIONS FOR USING THIS EXAMPLE
===================================

This directory contains an example run of an accretion outburst/quiescent cooling. 

1. This directory may be copied anywhere on your computer. Note, however, that you need write permission to dStar/data to run as the code caches datafiles there.  I'll change this in the near future to allow you to set the cache location (useful for running on multi-user machines).

2. Edit make/makefile.  Replace the following macros with the correct paths

        DSTAR_LIB_DIR = /path/to/dStar/lib
        DSTAR_INC = /path/to/dStar/include

3. Edit src/run.f.  Replace this variable with the correct path

        character(len=*), parameter :: my_dStar_dir = '/path/to/local/dStar'

4. Change the parameters as desired in the `inlist`.  Note that 

        maximum_end_time = 4320000.0
            ! seconds  (this is a 50-day outburst)

    and
            
            Mdot = 1.0e17
    
    refer to the accretion outburst.  To set the controls for the quiescent cooling, edit these lines in `run.f`:
    
        ! Now set Mdot = 0
        s% Mdot = 0.0_dp
        ! set the starting model profile to be one larger than the current model
        s% starting_number_for_profile = s% model + 1
        ! We can reset the start time to zero for convenience
        s% start_time = 0.0
        ! and we'll have a 1000 day quiescent period
        s% maximum_end_time = 8.64d8
        
    (Yes, I know this should be done better; it is a leftover from some experiments I was trying.  Will fix soon.)
    
5. Build the code: `./mk`
    
6. And run it: `./run_dStar`
    
Happy modeling!
