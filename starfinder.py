# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:41:37 2024

@author: Sam
"""
#import relevant modules
from astropy.table import Table
from astropy.io import fits 
import matplotlib.pyplot as plt
import matplotlib as mpl
import pylab as plt 
import numpy as np
import math
from scipy import stats

#Reading the table containing DR17 data
a = fits.open('dr17_dr3_McMillan_astroNN_rev1.fits')
hdu = a[1]
aid = a[1].data['APOGEE_ID_1']
ateff = a[1].data['TEFF']
alogg = a[1].data['LOGG']
asn = a[1].data['SNR']
arv = a[1].data['VHELIO_AVG']
afield = a[1].data['FIELD']
asflags = a[1].data['STARFLAGS']
aaflags = a[1].data['ASPCAPFLAGS']

#define working sample

mask = ( (asn >= 50) & (ateff > 3750) & (ateff < 5500) & (alogg > -1) & \
       (asflags != 'BAD' ) & (aaflags != 'BAD' ))

ids = aid[mask]
mlogg = alogg[mask]
marv = arv[mask]
mafield = afield[mask]

#find the Ursa Minor afield heading
#print(sorted(set(mafield)))

#Identify members of each dwarf galaxy

Carina = ((mlogg < 3.0) & (209.7 < marv) & (marv < 236.1) & (mafield=='CARINA'))
Draco = ((mlogg < 3.0) & (-309.2 < marv) & (marv < -272.8) & (mafield=='DRACO'))
Fornax = ((mlogg < 3.0) & (31.9 < marv) & (marv < 78.7) & (mafield=='FORNAX'))
Sculptor = ((mlogg < 3.0) & (93.0 < marv) & (marv < 129.8) & (mafield=='SCULPTOR'))
Sextans = ((mlogg < 3.0) & (208.4 < marv) & (marv < 240) & (mafield=='SEXTANS'))
Ursa_Minor = ((mlogg < 3.0) & (-265.9 < marv) & (marv < -227.9) & (mafield=='URMINOR'))

                   
CarinaIds = ids[Carina]
print('Carina Objects:',len(CarinaIds))
DracoIds = ids[Draco]
print('Draco Objects:',len(DracoIds))
FornaxIds = ids[Fornax]
print('Fornax Objects:',len(FornaxIds))
SculptorIds = ids[Sculptor]
print('Sculptor Objects:',len(SculptorIds))
SextansIds = ids[Sextans]
print('Sextans Objects:',len(SextansIds))
Ursa_MinorIds = ids[Ursa_Minor]
print('Ursa Minor Objects:',len(Ursa_MinorIds))

#combine all of these into 1 table

DwarfIds = np.concat((CarinaIds, DracoIds, FornaxIds, SculptorIds, SextansIds, Ursa_MinorIds))
print(len(DwarfIds))

#find the the ids of each of the dwarf galaxy stars

mask = np.isin(aid, DwarfIds)

#create new fits file to write Dwarf galaxies to

selected_data = hdu.data[mask]
image_hdu = fits.BinTableHDU(data=selected_data, name='Dwarf Galaxy Data')
print(image_hdu)
image_hdu.writeto('DwarfGalaxies.fits', overwrite=True)