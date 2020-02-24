INSTRUCTIONS FOR USING THIS EXAMPLE
===================================

This directory contains an example run of an accretion outburst/quiescent cooling with a more sophisticated input of the accretion history.

1.  This directory may be copied anywhere on your computer. Note, however, that you need write permission to dStar/data to run as the code caches datafiles there.  I'll change this in the near future to allow you to set the cache location (useful for running on multi-user machines).

2.  Edit make/makefile.  Replace the following macros with the correct paths

        DSTAR_LIB_DIR = /path/to/dStar/lib
        DSTAR_INC = /path/to/dStar/include

3.  Change the parameters as desired in the `inlist`.
    
4.  Build the code: `./mk`
    
5.  And run it: `./run_dStar -D<path to dStar> -I<inlist name>`
    The command now takes two optional arguments, the path to the root directory of dStar and the name of the inlist file.  You no longer need to edit the src code for that.  The default values are
    
        default_dStar_dir = '../../dStar'
        default_inlist_file = 'inlist'
        
6.  The accretion history is now contained in the file 'accretion_history'. As in Fortran, a '!' signifies that the rest of the line is to be ignored. The file consists of name/value pairs; if the name 'epochs' is encountered, the remaining lines are presumed to be a table of times and accretion rates.

    There are 4 values controlling the accretion history

    1.  number_cycles: allows for a repeating pattern

    2.  Mdot_scale: useful if the table is in Eddington values and need to 
        be converted to cgs scale

    3.  time_scale: for converting time to solar days

    4.  columns: indicates order of columns, default order is (time, Mdot)
    
    For n lines read and m cycles, we repeat lines 1:n-1 m times. The times are stored in `epoch_boundaries` at indices 0:(n-1)-1, n-1:2(n-1)-1, ..., (m-1)•(n-1):m(n-1)-1; if the table has times t1,...,tn, then with Delta = tn-t1 the stored times are t1, t2,...,t(n-1),t1+Delta,...,t(n-1)+Delta,t_1+2•Delta,... t(n-1)+(m-1)•Delta. The final time is stored at (m•(n-1)) = tn+(m-1)*Delta.

    The accretion rates are stored in `epoch_Mdots` at indices 1:n-1, (n-1)+1:2(n-1), ..., (m-1)*(n-1)+1:m(n-1).
    
Happy modeling!
