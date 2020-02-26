INSTRUCTIONS FOR USING THIS EXAMPLE
===================================

This directory contains an example run of an accretion outburst/quiescent cooling for which we wish to compare to an observed lightcurve (KS 1731-260; parameters are from Brown & Cumming [2009]). 

1. This directory may be copied anywhere on your computer. Note, however, that you need write permission to dStar/data to run as the code caches datafiles there.  I'll change this in the near future to allow you to set the cache location (useful for running on multi-user machines).

2. Edit make/makefile.  Replace the following macros with the correct paths

        DSTAR_LIB_DIR = /path/to/dStar/lib
        DSTAR_INC = /path/to/dStar/include

3. Change the parameters as desired in the `inlist`.  In particular,
    
        write_interval_for_terminal = 1000
        write_interval_for_history = 1000
        write_interval_for_profile = 1000
        starting_number_for_profile = 1

    These ensure that none of the logs are written to stdout, and that no profile and history files are made.  We want to monitor the effective temperature when there were observations.  The following settings do this; the end of the outburst is at t = 0.0 d.
    
        ! integration epochs
        number_epochs = 9
        basic_epoch_Mdots = 1.0e17,8*0.0
        basic_epoch_boundaries = -4383.0,0.0,65.1,235.7,751.6,929.5,1500.5,1570.4,1595.4,3039.7
        
    We modify `run.f` to compute a metric for how well the model fits the observed lightcurve.
    
        pred_Teff = s% Teff_monitor(2:)/1.0e6
        eV_to_MK = 1.602176565e-12_dp/boltzmann/1.0e6
        ! observed effective temperatures (eV) and uncertainties
        obs_Teff = [103.2,88.9,75.5,73.3,71.0,66.0,70.3,63.1] * eV_to_MK
        obs_Teff_sig = [1.7,1.3,2.2,2.3,1.8,4.5,2.1,2.1] * eV_to_MK
        chi2 = sum((pred_Teff-obs_Teff)**2 / obs_Teff_sig**2 )
        
    Look at the code for more details.
    
4. Build the code: `./mk`
    
5.  And run it: `./run_dStar -D<path to dStar> -I<inlist name>`
    The command now takes two optional arguments, the path to the root directory of dStar and the name of the inlist file.  You no longer need to edit the src code for that.  The default values are
    
        default_dStar_dir = '../../dStar'
        default_inlist_file = 'inlist'
    
Happy modeling!
