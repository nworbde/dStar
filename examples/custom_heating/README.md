INSTRUCTIONS FOR USING THIS EXAMPLE
===================================

This directory contains an example run of an accretion outburst/quiescent cooling with a custom routine to set Qimp and the heating.

1.  This directory may be copied anywhere on your computer. Note, however, that you need write permission to dStar/data to run as the code caches datafiles there.  I'll change this in the near future to allow you to set the cache location (useful for running on multi-user machines).

2.  Edit make/makefile.  Replace the following macros with the correct paths

        DSTAR_LIB_DIR = /path/to/dStar/lib
        DSTAR_INC = /path/to/dStar/include

3.  Change the parameters as desired in the `inlist`.  In particular, we want to monitor the observed effective temperature when there were observations. The following settings do this; the end of the outburst is at t = 0.0 d.
    
        ! integration epochs
        number_epochs = 12
        epoch_Mdots = 2.0e17,11*0.0
        epoch_boundaries = -912.0,0.0,30.0,56.0,100.0,180.0,300.0,560.0,1000.0,1800.0,3000.0,5600.0,10000.0

    To use the alternate routines, we turn them on in the inlist
   
        ! other routines
        use_other_set_Qimp = .TRUE.
        use_other_set_heating = .TRUE.
       
    We can pass extra controls to the routines using the `extra_*_controls`.

        ! extra controls for hook routines
        ! defined here
        ! 1. extra heating from pion -> neutrino
        ! 2. Q in the pasta
        ! 3.-4. density limits for extra heating
        extra_real_controls = 4.0,20.0,1.0e12,1.0e13
       
4.  We add a file `alt_micro.f` and put this in the makefile as well.  The interface for the routines is specified in `NScool/public/NScool_def.f`.

5. We add the module `alt_micro` to `run.f` and set the function pointers to our custom routines

        call get_NScool_info_ptr(NScool_id,s,ierr)
        s% other_set_Qimp => alt_Qimp
        s% other_set_heating => alt_heating


6.  Finally, we add in `run.f` some code to output the lightcurve
    
        pred_Teff = s% Teff_monitor(2:)/1.0e6    
        write(output_unit,*)
        write(output_unit,'(a7,a6)') 'time','Teff'
        write(output_unit,'(a7,a6)') '[d]','[MK]'
        write(output_unit,'(13("-"))')
        do i = 1, 11
            write(output_unit,'(f7.1,f6.3)')  &
            & s% t_monitor(i+1), pred_Teff(i)
        end do
        
    Look at the code for more details.
    
7.  Build the code: `./mk`
    
8.  And run it: `./run_dStar -D<path to dStar> -I<inlist name>`
    The command now takes two optional arguments, the path to the root directory of dStar and the name of the inlist file.  You no longer need to edit the src code for that.  The default values are
    
        default_dStar_dir = '../../dStar'
        default_inlist_file = 'inlist'
    

    
Happy modeling!
