SETTING THE CRUST COMPOSITION FROM AN ABUNTIME FORMAT FILE
==========================================================

By default, a composition table using the Haensel & Zdunik (1990) composition is created when the module `dStar_crust` is built. To install a composition in abuntime format, follow the following directions.

1. Place the abuntime file in the 'data' directory.
2. Run the `preprocess_composition` script. This script takes three arguments: the name of the abuntime file (preferably without extensions, as that name will be used as the "stem" for other files); a parameter `-d` that sets the increment in `lg(P)` between entries in the composition table; a parameter `-t` that sets the abundance (relative to the maximum at each point) for an isotope to be included in the composition table; and a parameter `-p` that sets a minimum abundance when printing out a summary of the table. Thus, for example,

        21:20:47$ ./preprocess_composition net_rp -d 0.005 -t 0.001 -p 0.01

    Here `net_rp` is the name of the abuntime file, the resulting preprocessed will have entries separated by roughly 0.005 in `lg(P)`, the table will be pruned to contain any isotopes with abundance exceeding 0.001 of the maximum at some location.
    
3.  There will be three files generated: if the abuntime file is `net_rp` as in this example, the files will be named `net_rp.bin`, `net_rp_isos` and `net_rp_summary`. The first of these is a binary format composition file used by `dStar`. The script copies it to the directory `dStar/data/crust_data/`. The second and third of these files contain a list of isotopes in the composition file and a summary of the composition as a function of pressure in the crust.

        21:43:37$ ./preprocess_composition net_rp -d 0.005 -t 0.01 -p 0.01
        nucchem_init: Loading nuclib from ../../data//nucchem
        nucchem_init: Retrieved 6342 nuclides. Writing nuclide dictionary
        process_abuntime: reading ../data/net_rp
        read_abuntime: nion =  1400
        ..............................
        read_abuntime: 717 zones accepted; 29879 zones rejected
        reduce_abuntime: setting min abundance to  1.00E-02
        extend_abuntime: extending lgP range from  26.835-- 30.598 to  22.000-- 33.498
        process_abuntime: writing cache to ../data/net_rp.bin
        nucchem_init: Loading nuclib from ../../data//nucchem
        nucchem_init: Retrieved 6342 nuclides. Writing nuclide dictionary
        print_composition: max number of isotopes with Y >  1.0000E-02*Ymax is   33 at lgP =  27.111
        
The file `net_rp_isos` looks like

           n
         o24
         f27
        ne28
        ne30
        ne32
        mg28
        mg30
        mg31
        mg32
        ...
         se76
        mg44 
        mg48 
        si50 
        si54 
        s52  
        s56  
        s60  
        ar66 
        ti88 
        
where you will note the last few isotopes are from Haensel and Zdunik (1990) and come from extending the composition to nuclear densities.  The summary table looks like

        22.000       n  0.000E+00  si28  1.749E-03  si30  7.076E-05   p31  9.091E-05   s32  3.463E-03  ...
         22.005       n  0.000E+00  si28  1.749E-03  si30  7.076E-05   p31  9.091E-05   s32  3.463E-03 ...
         22.010       n  0.000E+00  si28  1.749E-03  si30  7.076E-05   p31  9.091E-05   s32  3.463E-03 ...
         22.015       n  0.000E+00  si28  1.749E-03  si30  7.076E-05   p31  9.091E-05   s32  3.463E-03 ...
         22.020       n  0.000E+00  si28  1.749E-03  si30  7.076E-05   p31  9.091E-05   s32  3.463E-03 ...
         22.025       n  0.000E+00  si28  1.749E-03  si30  7.076E-05   p31  9.091E-05   s32  3.463E-03 ...
        ...
        33.468       n  8.036E-01 ti88   2.232E-03
        33.473       n  8.036E-01 ti88   2.232E-03
        33.478       n  8.036E-01 ti88   2.232E-03
        33.483       n  8.036E-01 ti88   2.232E-03
        33.488       n  8.036E-01 ti88   2.232E-03
        33.493       n  8.036E-01 ti88   2.232E-03
        33.498       n  8.036E-01 ti88   2.232E-03
        
where the first column is `lg(P)` followed by the most prominent isotopes at that pressure. To use the composition, specify the following variables in the inlist.

        real(dp) :: crust_reference_temperature ! default is 1.0e8K
        character(len=16) :: crust_composition  ! default is 'HZ90'
        
Here the string (16 character max) `crust_composition` is the name of the abuntime file.


