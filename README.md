#dStar

A collection of modules for computing neutron star structure and evolution.

##What's new

You can now specify a run with a number of accretion "epochs": distinct periods of time with a different accretion rate.  For example, suppose you wish to accrete at 1.5e17 g/s for 1000 d and then cool for 5000 d.  In the inlist, you would set the following flags.

      start_time = 0.0
      number_epochs = 2
      epoch_Mdots = 1.5e17, 0.0
      epoch_end_times = 1000.0, 6000.0

You can also use this to specify times at which you want the surface effective temperature recorded.  For example, suppose we want the surface effective temperature 50 d, 100 d, 500 d, 1000 d, 2000 d, and 5000 d after then end of the outburst in the above example.  We would then put the following in the inlist.

    start_time = 0.0
    number_epochs = 7
    epoch_Mdots = 1.5e17,6*0.0
    epoch_end_times = 1000.0, 1050.0, 1100.0, 1500.0, 2000.0, 3000.0, 6000.0

The structure pointer now contains arrays `t_monitor` and `Teff_monitor` that contain the epoch end times (in days) and the observer-frame effective temperature (in K) at those times.  This facilitates comparison with observations.

##Dependencies
  * [MESA](http://mesa.sourceforge.net): `dStar` makes use of the `MESA` numerical, utility, and equation of state libraries.
  * [MESA SDK](http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk): the compilation of both `MESA` and `dStar` has been tested using a specific build environment.

This version of `dStar` has been tested with `MESA` version 7503 and the 2014 December 12 (or later) version of the `MESA SDK`.

##How to install
  1. Follow the instructions on the MESA website to build a working version of `MESA`, and ensure that the environment variable `MESA_DIR` points to that directory.
  2. After forking or cloning the dStar repository, go into the top-level directory and type `./install`.  If you had compiled it previously, you should do a `./clean` first to force a clean install.

###What the installation does
For each module, the install script

  1. Downloads, verifies, and installs to the `data` directory the necessary data files
  2. Compiles each module as a library.
  3. Performs a test of each module and compares the output against a sample.  A deviation from allowed tolerances results in a install failure.
  4. Installs the library and module files into the top-level `install` and `lib` direcories.

##How to use
For each module, look in the `test` directory for an example of how to run the module. The primary module is `NScool`.

For a basic example of how to run a neutron star model over an accretion/quiescent cycle, copy `examples/basic_run` and follow the instruction in the `README.md` file in that directory.

##How to cite
If you do use `dStar`, we'd appreciate a citation! `dStar` is listed in the Astrophysics Source Code Library [ascl:1505.034](http://ascl.net/1505.034) and can be cited as, e.g., 
    
    Brown, E. F. 2015, Astrophysics Source Code Library, ascl:1505.034

A bibliographic entry can be obtained from [ADS](http://adsabs.harvard.edu/abs/2015ascl.soft05034B).


##Upcoming improvements
  1. add load/save options for models
  2. ability to generate an atmosphere model with an arbitrary composition.
  3. Add environment variable pointer to root directory and allow cache directories to be in a user-specified location.

  
