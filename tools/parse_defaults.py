from __future__ import print_function, division
from os import environ

default_file = environ['DSTAR_DIR']+'/NScool/defaults/controls.defaults'

defaults = {}

with open(default_file,'r') as f:
    for line in f:
        # strip comments and whitespace: everthing to right of ! is ignored
        toks = line.split('!')[0].split()
        if len(toks) == 0: continue
        # now smash toks back together and split on '='
        toks = ''.join(toks).split('=')
        k, v = toks[0],toks[1]
        defaults[k] = v

for k,v in defaults.iteritems():
    print(k+' = '+v)
