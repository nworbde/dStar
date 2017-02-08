INSTRUCTIONS FOR USING THIS EXAMPLE
===================================

This directory contains an example run of an accretion outburst/quiescent cooling with command-line options to set the shallow heating, and the core mass and radius.

1.  This directory may be copied anywhere on your computer. Note, however, that you need write permission to dStar/data to run as the code caches datafiles there.  I'll change this at some point to allow you to set the cache location (useful for running on multi-user machines).

2.  Edit make/makefile.  Replace the following macros with the correct paths

        DSTAR_LIB_DIR = /path/to/dStar/lib
        DSTAR_INC = /path/to/dStar/include

3.  Change the parameters as desired in the `inlist`.  In particular, we want to monitor the heating of the crust over the course of an outburst/quiescent cycle. The following settings do this; the end of the outburst is at t = 0.0 d.
    
        ! integration epochs
        number_epochs = 16
        epoch_Mdots = 10*6.7e16,6*0.0
        epoch_boundaries = -20.0,-18.0,-16.0,-14.0,-12.0,-10.0,-8.0,-6.0,-4.0,-2.0, 0.0,3.0,10.0,30.0,100.0,300.0,600.0
    
4.  Build the code: `./mk`
    
5.  And run it: `./run_dStar -D<path to dStar> -I<inlist name> -Q<shallow heating in Mev> -M<core mass in solar units> -R<core radius in km>`
    The command now takes five optional arguments, the path to the root directory of dStar and the name of the inlist file, as well as the amount of shallow heating and the core mass and radius. The default values are
    
        default_dStar_dir = '../../dStar'
        default_inlist_file = 'inlist'
        default_Q_heating_shallow = 0.0 ! MeV
        default_core_mass = 1.6 ! Msun
        default_core_radius = 11.0 ! km

Happy modeling!
