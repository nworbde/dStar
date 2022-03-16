import numpy as np

class dStarData:
    """
    Holds dStar data.  Each item is stored as a class attribute.
    The class is initialized with a dictionary of values.
    """
    def __init__(self,d):
        trans = str.maketrans('( ','__',').&=-+*%^#@!,/?<>;:~[]')
        for k, v in d.items():
            k = k.translate(trans)
            setattr(self,k,v)

class dStarReader:
    """
    Lightweight reader for dStar output files.
    
    Example
    =======
    
    In [1]: from reader import dStarReader

    In [2]: f = dStarReader('history.data')

    In [3]: f.headers
    Out[3]: ['gravity', 'core_mass', 'core_radius', 'core_temp']

    In [4]: f.header.gravity
    Out[4]: 202235100000000.0

    In [5]: f.columns
    Out[5]: ['lgLnuc', 'time', 'lgTeff', 'lgLsurf', 'lgLnu', 'lgMdot', 'model']

    In [6]: t = f.data.time

    In [7]: t
    Out[7]: 
    array([  0.00000000e+00,   1.15740700e-11,   4.62963000e-11,
             1.50463000e-10,   4.62963000e-10,   1.40046300e-09,
             4.21296300e-09,   1.26504600e-08,   3.79629600e-08,
             1.13900500e-07,   3.41713000e-07,   1.02515000e-06,
             3.07546300e-06,   9.22640000e-06,   2.76792100e-05,
             8.30376500e-05,   2.49113000e-04,   7.47338900e-04,
             2.24201700e-03,   6.72605000e-03,   2.01781500e-02,
             6.05344500e-02,   1.81603400e-01,   5.44810100e-01,
             1.63443000e+00,   4.90329100e+00,   1.47098700e+01,
             3.61640100e+01,   7.48279600e+01,   1.49495200e+02,
             2.94449100e+02,   5.78937500e+02,   1.12066100e+03,
             1.68763400e+03,   3.34658200e+03,   6.19151200e+03,
             1.00000000e+04,   0.00000000e+00,   1.00000000e-02,
             4.00000000e-02,   1.30000000e-01,   4.00000000e-01,
             1.21000000e+00,   3.64000000e+00,   1.09300000e+01,
             3.28000000e+01,   8.15802200e+01,   1.77517500e+02,
             3.62505400e+02,   7.10275100e+02,   1.35315800e+03,
             2.11398300e+03,   3.05794800e+03,   4.53322900e+03,
             7.30039400e+03,   1.00000000e+04])
    
    """
    
    def _set_type(self,k):
        return int if 'model' in k else float
    
    def __init__(self,file,meta_line=4,header_lines=8):
        # extra the meta data from the file
        with open(file) as f:
            lno = 0
            for line in f:
                lno+=1
                if lno == meta_line:
                    self._meta_head = line.split()
                    self._meta_data = np.fromstring(f.readline(),sep=' ')
                    break

        # stuff the meta data into a dictionary
        self._info = { k:v for k,v in zip(self._meta_head,self._meta_data)}
        # if one of the headers is 'model', keep this as an integer
        if 'model' in self._info:
            self._info['model'] = int(self._info['model'])
        self._header = dStarData(self._info)
        
        # now read in main section of data and convert to a recarry view
        self._arr = np.genfromtxt(file,skip_header=header_lines,names=True,dtype=None).view(np.recarray)
        
        # if this has a time series,
        # only keep models where t[i] - t[i-1] > 0
        if 'time' in self._arr.dtype.names:
            self._kept_models = [0]
            for i in range(1,self._arr.shape[0]):
                if self._arr.time[i] > self._arr.time[self._kept_models[-1]]:
                    self._kept_models.append(i)

    @property
    def header(self):
        """
        Metadata about the run.  To see the entries, use the headers attribute.
        """
        return self._header
    
    @property
    def data(self):
        """
        Data from the run, stored as a recarray.  To see the available entries, use the columns
        attribute.
        """
        return self._arr

    @property
    def headers(self):
        """
        List of available metadata
        """
        return self._meta_head
    
    @property
    def columns(self):
        """
        List of available column data
        """
        return self._arr.dtype.fields.keys()
    
    @property
    def num_lines(self):
        """
        number of lines in the file
        """
        return self._arr.shape[0]
    
    @property
    def increasing(self):
        return self._kept_models

class dStarCrustReader:
    """
    Lightweight reader for dStar profiles, i.e., crust properties as a function of position and time.
    
    Example
    =======
    
    In [1]: from reader import dStarCrustReader

    In [2]: f = dStarCrustReader('LOGS')

    In [3]: f.num_zones
    Out[3]: 253

    In [4]: f.num_models
    Out[4]: 37

    In [5]: f.parameters
    Out[5]: 
    ['eps_nu',
     'dm',
     'Xn',
     'temperature',
     'zone',
     'area',
     'luminosity',
     'K',
     'density',
     'gravity',
     'abar',
     'pressure',
     'ePhi',
     'mass',
     'eps_nuc',
     'zbar',
     'cp',
     'Qimp',
     'Gamma',
     'eLambda']

    In [6]: crust = f.crust

    In [7]: crust.temperature.shape
    Out[7]: (37, 253)

    In [8]: crust.temperature[:,0]  # this is the temperature in zone 0 (outer zone) for models 0...36.
    Out[8]: 
    array([  9.00452300e+07,   9.00452300e+07,   9.00452300e+07,
             9.00452300e+07,   9.00452300e+07,   9.00452300e+07,
             9.00452300e+07,   9.00452300e+07,   9.00452300e+07,
             9.00452300e+07,   9.00452300e+07,   9.00452300e+07,
             9.00452400e+07,   9.00452400e+07,   9.00452700e+07,
             9.00453500e+07,   9.00456000e+07,   9.00464000e+07,
             9.00491800e+07,   9.00597000e+07,   9.00992700e+07,
             9.02363600e+07,   9.06833300e+07,   9.20956400e+07,
             9.64892800e+07,   1.09370400e+08,   1.37259900e+08,
             1.69301000e+08,   1.97342000e+08,   2.23619000e+08,
             2.47974800e+08,   2.70005900e+08,   2.87756400e+08,
             2.96207400e+08,   3.04138900e+08,   3.05921900e+08,
             3.06066700e+08])   
    """

    def __init__(self,directory,**kwargs):
        # read the list of profiles
        profiles_list = directory+'/profiles'
        self._profiles = np.genfromtxt(profiles_list,skip_header=1,dtype=[('model','i8'), ('time','f8')]).view(np.recarray)
        self._models = self._profiles.model
        self._times = self._profiles.time

        dt = 0.0
        if 'dt' in kwargs:
            dt = kwargs.pop('dt')
            
        basename = 'profile'
        if 'basename' in kwargs:
            basename = kwargs.pop('basename')

        # only keep modes where t[i] - t[i-1] > dt
        kept_models = [0]
        for i in range(1,self.models.shape[0]):
            if self._times[i] > self._times[kept_models[-1]] + dt:
                kept_models.append(i)

        self._models = self._models[kept_models]
        self._times = self._times[kept_models]

        # set the shape
        self._nmodels = self._models.shape[0]
        # read the first profile to get layout
        d = dStarReader(directory+'/'+basename+str(self._models[0]))
        self._columns = d.columns
        self._nzones = d.num_lines
        self._crust = {}
        for col in d.columns:
            self._crust[col] = np.zeros((self._nmodels,self._nzones))
            self._crust[col][0,:] = d.data[col]

        for model in range(1,self._nmodels):
            d = dStarReader(directory+'/profile'+str(self._models[model]),**kwargs)
            for col in d.columns:
                self._crust[col][model,:] = d.data[col]
        
        self._crust = dStarData(self._crust)

    @property
    def parameters(self):
        """
        list of crust parameters
        """
        return self._columns

    @property
    def times(self):
        """
        times for each crust profile
        """
        return self._times
    
    @property
    def models(self):
        """
        model numbers for each stored crust profile
        """
        return self._models

    @property
    def crust(self):
        """
        returns a dStarData object contains the columns for the crust
        """
        return self._crust

    @property
    def num_models(self):
        """
        number of models stored
        """
        return self._nmodels
    
    @property
    def num_zones(self):
        """
        number of zones for each profile
        """
        return self._nzones
