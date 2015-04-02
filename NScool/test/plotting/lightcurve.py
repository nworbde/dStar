# script for dStar lightcurves

import numpy as np
from matplotlib import pyplot as plt

#dStar data
model, dstar_time, logmdot, dstar_logteff, log_lsurf, log_lnu, log_lnuc = np.loadtxt('../LOGS/history.data',unpack=True,skiprows=13)

dstar_teff = 10.0**(dstar_logteff)

time_no_mdot = []
teff_no_mdot = []
time_no_mdot_days = []
time_when_acc_halts = 0.

# redshift = (1-2GM/(rc^2))^(-1/2)
redshift = 1.0

# save times and teffs when there is no accretion in new arrays
for ii in range(len(logmdot)):
    if logmdot[ii] == 0 and time_when_acc_halts == 0:
        time_when_acc_halts = dstar_time[ii-1]
    if logmdot[ii] == 0:
        time_no_mdot.append(dstar_time[ii]-time_when_acc_halts)
        teff_no_mdot.append(dstar_teff[ii]*8.617E-5/redshift)

# convert time into log(t/d)
for jj in range(len(time_no_mdot)):
    time_no_mdot_days.append(time_no_mdot[jj]/86400.0)

fig = plt.figure()
ax = fig.add_subplot(1,1,1, adjustable='box')
dx = plt.gcf().set_size_inches(10,10)
ax.plot(time_no_mdot_days,teff_no_mdot,linewidth=2.0,color='k')

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2.0)

ax.set_xlabel('$\mathrm{log}(t) \ [\mathrm{days}]$',fontsize=30)
ax.set_ylabel('$T_{\mathrm{eff}}^{\infty} \ [\mathrm{eV}]$', fontsize=30)
ax.set_ylim(0,350)
ax.set_xscale('log',nonposx='clip')
ax.set_xlim(1,10000)
ax.tick_params(axis='x',length=8,width=2,which='major')
ax.tick_params(axis='y',length=8,width=2,which='major')
ax.tick_params(axis='x',length=5,width=2,which='minor')

plt.savefig('plots_out/lightcurve.pdf',format='pdf')
