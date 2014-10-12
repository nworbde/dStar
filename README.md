#dStar

A collection of modules for computing neutron star structure and evolution.

##Dependencies
  * [MESA](http://mesa.sourceforge.net): `dStar` makes use of the `MESA` numerical, utility, and equation of state libraries
  * [MESA SDK](http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk): the compilation of both `MESA` and `dStar` has been tested using a specific build environment.

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
**Under construction**
For each module, look in the `test` directory for an example of how to run the module.

##Upcoming improvements
  1. add load/save options for models
  