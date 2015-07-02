from __future__ import print_function, division
import numpy as np

class dStarReader:
    """
    Lightweight reader for dStar output files.
    
    Example
    =======
    
    In [1]: from reader import dStarReader

    In [2]: f = dStarReader('history.data')

    In [3]: f.list_headers()
    Out[3]: ['gravity', 'core_mass', 'core_radius', 'core_temp.']

    In [4]: f.get_header('gravity')
    Out[4]: 204558200000000.0

    In [5]: f.list_columns()
    Out[5]: ['lgLnuc', 'times', 'lgTeff', 'lgLsurf', 'lgLnu', 'lgMdot', 'model']

    In [6]: t = f.get_column('times')

    In [7]: t
    Out[7]: 
    array([  0.00000000e+00,   1.00000000e-06,   4.00000000e-06,
             1.30000000e-05,   4.00000000e-05,   1.21000000e-04,
             3.64000000e-04,   1.09300000e-03,   3.28000000e-03,
             9.84100000e-03,   2.95240000e-02,   8.85730000e-02,
             2.65720000e-01,   7.97161000e-01,   2.39148400e+00,
             7.17445300e+00,   2.15233600e+01,   6.45700800e+01,
             ...
             8.64000000e+08])
    
    """
    
    def __init__(self,file,meta_line=4,header_lines=8):
        # extra the meta data from the file
        with open(file) as f:
            lno = 0
            for line in f:
                lno+=1
                if lno == meta_line:
                    self._meta_head = line.split()
                    self._meta_data = np.fromstring(f.next(),sep=' ')
                    break

        # stuff the meta data into a dictionary
        self.info = { k:v for k,v in zip(self._meta_head,self._meta_data)}
        # if one of the headers is 'model', keep this as an integer
        if 'model' in self.info:
            self.info['model'] = int(self.info['model'])

        # now read in main section of data
        self._arr = np.genfromtxt(file,skip_header=header_lines,names=True,dtype=None)

    def get_header(self,name):
        return self.info[name]
    
    def get_column(self,name):
        return self._arr[:][name]

    def list_headers(self):
        return self._meta_head
    
    def list_columns(self):
        return self._arr.dtype.fields.keys()
        
    def num_lines(self):
        return self._arr.shape[0]

class dStarCrustReader:
    
    def __init__(self,directory,**kwargs):
        # read the list of profiles
        profiles_list = directory+'/profiles'
        self._profiles = np.genfromtxt(profiles_list,skip_header=1,dtype=[('model','i8'), ('time','f8')])
        self._models = self._profiles['model']
        self._times = self._profiles['time']
        
        dt = 0.0
        if 'dt' in kwargs:
            dt = kwargs.pop('dt')

        # only keep modes where t[i] - t[i-1] > dt
        kept_models = [0]
        for i in range(1,self._models.shape[0]):
            if self._times[i] > self._times[kept_models[-1]] + dt:
                kept_models.append(i)

        self._models = self._models[kept_models]
        self._times = self._times[kept_models]
        # read the first profile to get layout
        nmodels = self._models.shape[0]
        d = dStarReader(directory+'/profile'+str(self._models[0]))
        nzones = d.num_lines()
        self._crust = {}
        for col in d.list_columns():
            self._crust[col] = np.zeros((nmodels,nzones))
            self._crust[col][0,:] = d.get_column(col)
        
        for model in range(1,self._models.shape[0]):
            d = dStarReader(directory+'/profile'+str(self._models[model]),**kwargs)
            for col in d.list_columns():
                self._crust[col][model,:] = d.get_column(col)

    def parameters(self):
        return self._crust.keys()
    
    def times(self):
        return self._times
    
    def models(self):
        return self._models
    
    def get_crust_parameter(self,name):
        return self._crust[name]
