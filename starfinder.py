# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:41:37 2024

@author: Sam
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

#Reading the table containing DR17 data
a = fits.open('dr17_dr3_McMillan_astroNN_rev1.fits')
aid = a[1].data['APOGEE_ID_1']
agid = a[1].data['GAIAEDR3_SOURCE_ID']
ara = a[1].data['RA']
adec = a[1].data['DEC']
aglon = a[1].data['GLON']
aglat = a[1].data['GLAT']
avscat = a[1].data['VSCATTER']
anvis = a[1].data['NVISITS']
ateff = a[1].data['TEFF']
ateff_e = a[1].data['TEFF_ERR']
alogg = a[1].data['LOGG']
alogg_e = a[1].data['LOGG_ERR']
asn = a[1].data['SNR']
arv = a[1].data['VHELIO_AVG']
afield = a[1].data['FIELD']

ajmag = a[1].data['J']
ahmag = a[1].data['H']
akmag = a[1].data['K']

afe = a[1].data['Fe_H']
amg = a[1].data['Mg_Fe']
ac = a[1].data['C_Fe']
ao = a[1].data['O_Fe']
an = a[1].data['N_Fe']
aca = a[1].data['Ca_Fe']
asi = a[1].data['Si_Fe']
ani = a[1].data['Ni_Fe']
aal = a[1].data['Al_Fe']
ati = a[1].data['Ti_Fe']
aco = a[1].data['Co_Fe']
a_s = a[1].data['S_Fe']
ak = a[1].data['K_Fe']
ap = a[1].data['P_Fe']
ace = a[1].data['Ce_Fe']
av = a[1].data['V_Fe']
acr = a[1].data['Cr_Fe']
ana = a[1].data['Na_Fe']
amn = a[1].data['Mn_Fe']
acu = a[1].data['Cu_Fe']
asflag = a[1].data['STARFLAG']
aaflag = a[1].data['ASPCAPFLAG']
asflags = a[1].data['STARFLAGS']
aaflags = a[1].data['ASPCAPFLAGS']
aeal = a[1].data['Al_FE_ERR']
aemg = a[1].data['MG_FE_ERR']
aemn = a[1].data['MN_FE_ERR']

ajmag = a[1].data['J']
ahmag = a[1].data['H']
akmag = a[1].data['K']
aa_k = a[1].data['AK_TARG']
aa_j = 2.5*aa_k
#ai36 = a[1].data['IRAC_3_6']
#ai45 = a[1].data['IRAC_4_5']
#ai80 = a[1].data['IRAC_8_0']

adist = a[1].data['weighted_dist']
aedist = a[1].data['weighted_dist_error']
ae = a[1].data['ecc']
#ae_err = a[1].data['e_err']
azmax = a[1].data['zmax']
#azmax_err = a[1].data['zmax_err']
arperi = a[1].data['rperi']
#arperi_err = a[1].data['rperi_err']
ajr = a[1].data['jr']
#ajr_err = a[1].data['jr_err']
aLz = a[1].data['Lz']
#aLz_err = a[1].data['Lz_err']
ajz = a[1].data['jz']
#ajz_err = a[1].data['jz_err']
aEnergy = a[1].data['energy']
#aEnergy_err = a[1].data['Energy_err']


print(len(aid))

#Calculate X, Y, Z, and Rgc

xsun = 8.e3

axx = (xsun - adist*np.cos(aglon*np.pi/180.)*np.cos(aglat*np.pi/180.))/1.e3
ayy = adist*np.sin(aglon*np.pi/180.)*np.cos(aglat*np.pi/180.)/1.e3
azz = adist*np.sin(aglat*np.pi/180.)/1.e3
aRgc = np.sqrt(axx**2. + ayy**2. + azz**2.)

print(len(afe))

#define working sample

mask = ( (asn >= 50) & (ateff > 3750) & (ateff < 5500) & (alogg > -1) & \
       (asflags != 'BAD' ) & (aaflags != 'BAD' ))

ids = aid[mask]
mlogg = alogg[mask]
marv = arv[mask]
mafield = afield[mask]
print(len(ids))

#Identify members of each dwarf galaxy

Carina = ((mlogg < 3.0) & (209.7 < marv) & (marv < 236.1) & (mafield=='CARINA'))
Draco = ((mlogg < 3.0) & (-309.2 < marv) & (marv < -272.8) & (mafield=='DRACO'))
Fornax = ((mlogg < 3.0) & (31.9 < marv) & (marv < 78.7) & (mafield=='FORNAX'))
Sculptor = ((mlogg < 3.0) & (93.0 < marv) & (marv < 129.8) & (mafield=='SCULPTOR'))
Sextans = ((mlogg < 3.0) & (208.4 < marv) & (marv < 240) & (mafield=='SEXTANS'))
Ursa_Minor = ((mlogg < 3.0) & (-265.9 < marv) & (marv < -227.9) & (mafield=='URSAMINOR'))

                   
CarinaIds = ids[Carina]
print(len(CarinaIds))
DracoIds = ids[Draco]
print(len(DracoIds))
FornaxIds = ids[Fornax]
print(len(FornaxIds))
SculptorIds = ids[Sculptor]
print(len(SculptorIds))
SextansIds = ids[Sextans]
print(len(SextansIds))
Ursa_MinorIds = ids[Ursa_Minor]
print(len(Ursa_MinorIds))

