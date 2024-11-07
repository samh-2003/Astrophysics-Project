# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 11:58:01 2024

@author: samhe
"""

#READ THE DATA 
#
from astropy.table import Table
from astropy.io import fits 
import matplotlib.pyplot as plt
import matplotlib as mpl
import pylab as plt
import numpy as np
import math
from scipy import stats


# Get current size
fig_size = plt.rcParams["figure.figsize"]

# Set figure width to 16 and height to 12
fig_size[0] = 20
fig_size[1] = 20

#Figure settings (font sizes of ticks and labels)
labs=35
tcks=30

#Reading the table of Dwarf Galaxy Objects
a = fits.open('DwarfGalaxies.fits')
id = a[1].data['APOGEE_ID_1']
gid = a[1].data['GAIAEDR3_SOURCE_ID']
ra = a[1].data['RA']
dec = a[1].data['DEC']
glon = a[1].data['GLON']
glat = a[1].data['GLAT']
vscat = a[1].data['VSCATTER']
nvis = a[1].data['NVISITS']
teff = a[1].data['TEFF']
teff_e = a[1].data['TEFF_ERR']
logg = a[1].data['LOGG']
logg_e = a[1].data['LOGG_ERR']
sn = a[1].data['SNR']
rv = a[1].data['VHELIO_AVG']
field = a[1].data['FIELD']

jmag = a[1].data['J']
hmag = a[1].data['H']
kmag = a[1].data['K']

feh = a[1].data['Fe_H']
mgfe = a[1].data['Mg_Fe']
cfe = a[1].data['C_Fe']
ofe = a[1].data['O_Fe']
nfe = a[1].data['N_Fe']
cafe = a[1].data['Ca_Fe']
sife = a[1].data['Si_Fe']
nife = a[1].data['Ni_Fe']
alfe = a[1].data['Al_Fe']
tife = a[1].data['Ti_Fe']
cofe = a[1].data['Co_Fe']
sfe = a[1].data['S_Fe']
kfe = a[1].data['K_Fe']
pfe = a[1].data['P_Fe']
cefe = a[1].data['Ce_Fe']
vfe = a[1].data['V_Fe']
crfe = a[1].data['Cr_Fe']
nafe = a[1].data['Na_Fe']
mnfe = a[1].data['Mn_Fe']
cufe = a[1].data['Cu_Fe']
alfe_err = a[1].data['Al_FE_ERR']
aemg = a[1].data['MG_FE_ERR']
aemn = a[1].data['MN_FE_ERR']
mgmn = mgfe - mnfe
simn = sife - mnfe
mgmn_err = np.sqrt(aemg**2. + aemn)

jmag = a[1].data['J']
hmag = a[1].data['H']
kmag = a[1].data['K']
a_k = a[1].data['AK_TARG']
a_j = 2.5*a_k
#ai36 = a[1].data['IRAC_3_6']
#ai45 = a[1].data['IRAC_4_5']
#ai80 = a[1].data['IRAC_8_0']

adist = a[1].data['weighted_dist']
aedist = a[1].data['weighted_dist_error']
dist = adist/1.e3
dist_err = aedist/1.e3
ecc = a[1].data['ecc']
#ecc_err = a[1].data['e_err']
zmax = a[1].data['zmax']
#zmax_err = a[1].data['zmax_err']
rperi = a[1].data['rperi']
#rperi_err = a[1].data['rperi_err']
jr = a[1].data['jr']/1.e2*8.*220
#jr_err = a[1].data['jr_err']
lz = a[1].data['Lz']/1.e2*8.*220
#Lz_err = a[1].data['Lz_err']
jz = a[1].data['jz']/1.e2*8.*220
#jz_err = a[1].data['jz_err']
ener = a[1].data['energy']/1.e5
#ener_err = a[1].data['Energy_err']

#Calculate Galactocentric distances X, Y, Z, and Rgc

xsun = 8.e3

xx = (xsun - adist*np.cos(glon*np.pi/180.)*np.cos(glat*np.pi/180.))/1.e3
yy = adist*np.sin(glon*np.pi/180.)*np.cos(glat*np.pi/180.)/1.e3
zz = adist*np.sin(glat*np.pi/180.)/1.e3
Rgc = np.sqrt(xx**2. + yy**2. + zz**2.)

#Graphs
#Nitrogen Against Iron

fig_size[0] = 10

plt.scatter(feh, nfe, c='Black', alpha=0.75, s=30)
plt.xlabel('FE', size=labs)
plt.xticks(np.arange(-2.5, 0, step=0.5), fontsize=tcks)
plt.ylabel('N', size=labs)
plt.yticks(np.arange(-0.5, 2, step=0.5), fontsize=tcks)
#plt.legend(loc='lower right')
plt.show()




