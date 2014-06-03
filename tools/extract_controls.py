#!/usr/bin/env python

import re
import argparse

parser = argparse.ArgumentParser(description="reads a list of variable declarations and strips off everything to the left of the '::'.  This is useful for generating default listings and namelists")

parser.add_argument('-d','--default-assignment',help='append equals sign to variable name, followed by variable name with default_ prepended. Used for assigning default values',default=False,action='store_true')
parser.add_argument('-a','--add-ampersand',help='append line continuation character',action='store_true',default=False)
parser.add_argument('-s','--add-store',help='write in form of an assign to a member of storage class',action='store_true',default=False)
parser.add_argument('-n','--namelist',nargs=1,help='construct a namelist')
parser.add_argument('input_file',type=str,help='file containing variable definitions')

var = re.compile(r'::\s+(\w+)')

args = parser.parse_args()
eq = ''
amp = ''
nml = ''

if (args.default_assignment): eq = ' = '
if (args.add_ampersand): amp = ', &'
if (args.namelist): nml = args.namelist

controls_file = args.input_file
varlist=[]
with open(controls_file,'r') as f:
    for line in f:
        m = var.search(line)
        if m:
            varlist.append(m.group(1))

if args.default_assignment:
    for var in varlist:
        print '    {0} = default_{0}'.format(var)

if args.add_store:
    for var in varlist:
        print '    s% {0} = {0}'.format(var)

if nml:
    print 'namelist /{0}/ &'.format(nml[0])
    for var in varlist[0:-1]:
        print '    {0}, & '.format(var)
    print '    {0}'.format(varlist[-1])

if args.add_ampersand:
    for var in varlist:
        print '    {0}{1}{2}'.format(var,eq,amp)
