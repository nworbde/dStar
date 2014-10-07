#!/usr/bin/env python

"""
    parses superfluid data directory and write report of available options
"""

from os import listdir

proton_singlet = {}
neutron_singlet = {}
neutron_triplet = {}

known_sf = {}
known_sf['p1'] = proton_singlet
known_sf['n1'] = neutron_singlet
known_sf['n3'] = neutron_triplet


def parse_one(file):
    info = None
    with open(file,'r') as f:
        info = f.readline().rstrip()
        return info

datafiles = listdir('.')
headers = {'p1':'proton 1S0', 'n1':'neutron 1S0', 'n3':'neutron 3P2'}

for d in datafiles:
    base = d.split('.')[0]
    attr, sftype = base.split('_')
    info = parse_one(d)
    known_sf[sftype][attr] = info

print '!   {0:>8}    {1:<}'.format('attr','reference notes')

for sftype in ['p1','n1','n3']:
    print '\n!   {0:^16}\n!   ================'.format(headers[sftype])
    for entry in iter(known_sf[sftype]):
            print '!   {0:>8}    {1:<}'.format(entry,known_sf[sftype][entry])
