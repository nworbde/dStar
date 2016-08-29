from __future__ import print_function, division
from os import environ
from sys import stdout

class dStarRun:
    
    _num_arr_controls = 32
    
    def __init__(self,dStar_dir=''):
        self._controls = {}
        self._defaults = {}
        self._names = []
        if dStar_dir == '':
            try:
                dStar_dir = environ['DSTAR_DIR']
            except KeyError:
                print('supply the root directory of dStar, either as an argument or through the environment variable DSTAR_DIR')
                return

        _defaults_file = dStar_dir+'/NScool/defaults/controls.defaults'

        with open(_defaults_file,'r') as f:
            for line in f:
                # strip comments and whitespace: everthing to right of ! is ignored
                toks = line.split('!')[0].split()
                if len(toks) == 0: continue
                # now smash toks back together and split on '='
                toks = ''.join(toks).split('=')
                # extra care needed for arrays
                k, v = toks[0],toks[1]
                if '(:)' not in k:
                    self._defaults[k] = v
                    self._names.append(k)
                else:
                    self._names.append(k)
                    self._defaults[k] = str(dStarRun._num_arr_controls)+'*'+v
                    for i in range(1,dStarRun._num_arr_controls+1):
                        arr_k = k.replace('(:)','('+str(i)+')')
                        self._defaults[arr_k] = v

    def control_setting(self,key):
        v = None
        try:
            v = self._controls[key]
        except:
            try:
                v = self._defaults[key]
            except:
                print('no control named '+key)
        return v

    def print_namelist(self,filename=None,defaults=False):
        if filename is None:
            f = stdout
        else:
            try:
                f = open(filename,'w')
            except:
                print('unable to open '+filename)
                f = stdout
        f.write('&controls\n')
        if defaults:
            for control in self._names:
                f.write('    {0} = {1}\n'.format(control,self._defaults[control]))
        else:
            for control,value in self._controls.iteritems():
                f.write('    {0} = {1}\n'.format(control,value))
        f.write('/\n')

    def set_control(self,control,value):
        if control in self._defaults:
            self._controls[control] = value
        else:
            print('bad control value {0}'.format(control))
            raise KeyError(control)
        return self

