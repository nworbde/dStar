from __future__ import print_function, division
from os import environ

class dStarRun:
    
    def __init__(self,dStar_dir=''):
        self._controls = {}
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
                k, v = toks[0],toks[1]
                self._controls[k] = v

    def print_namelist(self):
        print('&controls')
        for control,value in self._controls.iteritems():
            print('    {0} = {1}'.format(control,value))
        print('/')

    def write_namelist_to_file(self,filename):
        with open(filename,'w') as f:
            f.write('&controls\n')
            for control,value in self._controls.iteritems():
                f.write('    {0} = {1}\n'.format(control,value))
            f.write('/\n')
    
    def set_control(self,control,value):
        if control in self._controls:
            self._controls[control] = value
        else:
            print('bad control value {0}'.format(control))
            raise KeyError(control)
