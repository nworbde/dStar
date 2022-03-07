# dStar

A collection of modules for computing neutron star structure and evolution.

## What's new

### Crust composition from `abuntime`-format files
As part of work with Rahul Jain and Hendrik Schatz, the crust composition is now specified by loading tables. Instructions are in the README file in the module `dStar_crust`. There are now scripts to preprocess files in an `abuntime` format (output from reaction network runs) into locally cached files of the composition.

### Neutron thermal conductivity
Alex Deibel wrote routines to compute the thermal conductivity of normal neutrons (Deibel et al. [2017], Astrophys. Jour. 839: 95). In addition, the thermal conductivity via phonons in the superfluid phase (Aguilera et al. [2009], Phys. Rev. Lett. 102: 091101) is computed. The neutron conductivity is not included by default, but can be activated by setting 

    use_sf_conductivity = .TRUE.
    use_nph_conductivity = .TRUE.
    use_nQ_conductivity = .TRUE.

in the inlist.

### Multiple accretion epochs
You can specify a run with a number of accretion "epochs": distinct periods of time with a different accretion rate.  For example, suppose you wish to accrete at 1.5e17 g/s for 1000 d (starting at t = 0 d) and then cool for 5000 d (that is, from t = 1000 d to t = 6000 d).  In the inlist, you would set the following flags.

      number_epochs = 2
      basic_epoch_Mdots = 1.5e17,0.0
      basic_epoch_boundaries = 0.0,1000.0,6000.0

You can also use this to specify times at which you want the surface effective temperature recorded.  For example, suppose we want the surface effective temperature 50 d, 100 d, 500 d, 1000 d, 2000 d, and 5000 d after the end of the outburst in the above example.  We would then put the following in the inlist.

    number_epochs = 7
    basic_epoch_Mdots = 1.5e17,6*0.0
    basic_epoch_boundaries = 0.0,1000.0,1050.0,1100.0,1500.0,2000.0,3000.0,6000.0

We can make this even more convenient by setting the end of the outburst at t = 0 d.

    number_epochs = 7
    basic_epoch_Mdots = 1.5e17,6*0.0
    basic_epoch_boundaries = -1000.0,0.0,50.0,100.0,500.0,1000.0,2000.0,5000.0

The number of "basic" epochs that can be set in the namelist is capped at 64. More complicated histories can be loaded from a file:

    load_epochs = .TRUE.
    epoch_datafile = 'accretion_history'

The `epoch_datafile`, in this example "accretion_history" has a header followed by a table of times and accretion rates. Within the header, `!` are used to indicate comments, and there are 4 control parameters, as indicated in the following example.

    number_cycles 2
    Mdot_scale 1.0e18  ! [g/s]
    time_scale 1.0     ! [d]
    columns 'time [d]'  'Mdot [Eddington]'

To repeat the outburst/quiescent cycle multiple times, set `number_cycles` > 1. The scale factors, `Mdot_scale` and `time_scale`, are multipliers for the table that convert to `g/s` and `d` for the accretion rate and time. Finally, the `columns` entry indicates the column order. Only the first 4 characters are examined to see if the column contains 'time' or 'Mdot'

The header ends with a line
    
    epochs

Subsequent lines are assumed to contain two columns of time and accretion rate. These columns are read into a temporary buffer, scaled by `Mdot_scale` and `time_scale`, repeated `number_cycles` times, and then used to integrate in time.  See the example `accretion` for an example.

### A tutorial on crust cooling
Look in the `examples/INT-16-2b-demo` directory for a demonstration of using this code that was presented in a talk given at the INT workshop 16-2b, "Phases of Dense Matter".  See the `README.md` file in that directory for instructions.

### Total mass, radius are now reported
Both the `history.data` and `profile` datafiles now list the total mass (solar units) and radius (kilometers) in the header.  By total, I mean the value at the top of the domain, not the photosphere.

### Time, observed effective temperature are stored
The structure pointer now contains arrays `t_monitor` and `Teff_monitor` that contain the epoch end times (in days) and the observer-frame effective temperature (in K) at the end of each epoch.  This facilitates comparison with observations.

### Command line options
Check out the examples `basic_run`, `fit_lightcurve`, and `custom_Tc_Qimp` to see how to add command line options to your run.

### Tools for reading output
Check out `tools/reader.py`. This contains a python class for reading the output history and profile data files. The quantities are stored as class members for easy access when analyzing results.

## Dependencies
  * [git-lfs](https://git-lfs.github.com): the data files are stored using `git lfs`, which needs to be installed prior to cloning the repository.
  * [MESA](http://mesa.sourceforge.net): `dStar` makes use of the `MESA` numerical, utility, and equation of state libraries. Because of changes to how `MESA` handles physical constants, **you must use version 15140 or later**.
  * [MESA SDK](http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk): the compilation of both `MESA` and `dStar` has been tested using a specific build environment.

This version of `dStar` has been tested with `MESA` version 15140 and the 2020 December 18 version of the `MESA SDK` for MacOS, `mesasdk-x86_64-macos-20.12.2`.

## How to install
  1. Follow the instructions on the MESA website to build a working version of `MESA`, and ensure that the environment variable `MESA_DIR` points to that directory.
  2. Install git lfs from https://git-lfs.github.com/ appropriate for the operating system you are using. Run the command: $ git lfs install. Reclone the dStar repository if necessary.
  3. After forking or cloning the dStar repository, go into the top-level directory and type `./install`.  If you had compiled it previously, you should do a `./clean` first to force a clean install.

### What the installation does
For each module, the install script

  1. Downloads, verifies, and installs to the `data` directory the necessary data files.  Note that the `md5` utility is used first for checking data integrity. This may be an issue on linux systems which prefer `md5sum`. If so, edit the `fetch_data` scripts.
  2. Compiles each module as a library.
  3. Performs a test of each module and compares the output against a sample.  A deviation from allowed tolerances results in a install failure.
  4. Installs the library and module files into the top-level `install` and `lib` directories.

## How to use
For each module, look in the `test` directory for an example of how to run the module. The primary module is `NScool`.

For a basic example of how to run a neutron star model over an accretion/quiescent cycle, copy `examples/basic_run` and follow the instruction in the `README.md` file in that directory.

## How to cite
If you do use `dStar`, we'd appreciate a citation! `dStar` is listed in the Astrophysics Source Code Library [ascl:1505.034](http://ascl.net/1505.034) and can be cited as, e.g., 
    
    Brown, E. F. 2015, Astrophysics Source Code Library, ascl:1505.034

A bibliographic entry can be obtained from [ADS](http://adsabs.harvard.edu/abs/2015ascl.soft05034B).

## Planned improvements
  1. add load/save options for models
  2. ability to generate an atmosphere model with an arbitrary composition.
  3. Add environment variable pointer to root directory and allow cache directories to be in a user-specified location.

  
