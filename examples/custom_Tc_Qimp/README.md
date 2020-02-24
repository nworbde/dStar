INSTRUCTIONS FOR USING THIS EXAMPLE
===================================

This directory contains an example run of an accretion outburst/quiescent cooling with a custom routine to set Qimp and the superfluid transition temperature.

1.  This directory may be copied anywhere on your computer. Note, however, that you need write permission to dStar/data to run as the code caches datafiles there.  I'll change this in the near future to allow you to set the cache location (useful for running on multi-user machines).

2.  Edit make/makefile.  Replace the following macros with the correct paths

        DSTAR_LIB_DIR = /path/to/dStar/lib
        DSTAR_INC = /path/to/dStar/include

3.  Change the parameters as desired in the `inlist`.  In particular, we want to monitor the observed effective temperature when there were observations. The following settings do this; the end of the outburst is at t = 0.0 d.
    
        ! integration epochs
        number_epochs = 9
        basic_epoch_Mdots = 1.0e17,8*0.0
        basic_epoch_boundaries = -4383.0,0.0,65.1,235.7,751.6,929.5,1500.5,1570.4,1595.4,3039.7

    To use the alternate routines, we turn them on in the inlist
   
        ! other routines
        use_other_set_Qimp = .TRUE.
        use_other_sf_critical_temperatures = .TRUE.
       
    We can pass extra controls to the routines using the `extra_*_controls`.

        ! extra controls for hook routines
        ! defined here
        ! 1. density at which to switch to high Q and neutron density at which to 
        !    switch to high Tc
        ! 2. Q in the pasta
        extra_real_controls = 4.0e13, 30.0
       
4.  We add a file `alt_micro.f` and put this in the makefile as well.  The interface for the routines is specified in `NScool/public/NScool_def.f`.

5. We add the module `alt_micro` to `run.f` and set the function pointers to our custom routines

        call get_NScool_info_ptr(NScool_id,s,ierr)
        s% other_set_Qimp => alt_Qimp
        s% other_sf_get_results => alt_sf


6.  Finally, we add in `run.f` some code to compute a metric for how well the model fits the observed lightcurve of KS1731-260.
    
        pred_Teff = s% Teff_monitor(2:)/1.0e6
        eV_to_MK = 1.602176565e-12_dp/boltzmann/1.0e6
        ! observed effective temperatures (eV) and uncertainties
        obs_Teff = [103.2,88.9,75.5,73.3,71.0,66.0,70.3,63.1] * eV_to_MK
        obs_Teff_sig = [1.7,1.3,2.2,2.3,1.8,4.5,2.1,2.1] * eV_to_MK
        chi2 = sum((pred_Teff-obs_Teff)**2 / obs_Teff_sig**2 )
        
    Look at the code for more details.
    
7.  Build the code: `./mk`
    
8.  And run it: `./run_dStar -D<path to dStar> -I<inlist name>`
    The command now takes two optional arguments, the path to the root directory of dStar and the name of the inlist file.  You no longer need to edit the src code for that.  The default values are
    
        default_dStar_dir = '../../dStar'
        default_inlist_file = 'inlist'
    

    
Happy modeling!
