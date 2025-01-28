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
labs=40
tcks=35

#Reading the table of Dwarf Galaxy Objects
a = fits.open('DwarfGalaxies.fits')
glon = a[1].data['GLON']
glat = a[1].data['GLAT']
field = a[1].data['FIELD']
feh = a[1].data['Fe_H']
mgfe = a[1].data['Mg_Fe']
cfe = a[1].data['C_Fe']
nfe = a[1].data['N_Fe']
alfe = a[1].data['Al_Fe']
adist = a[1].data['weighted_dist']

#Calculate Galactocentric distances X, Y, Z, and Rgc

xsun = 8.e3

xx = (xsun - adist*np.cos(glon*np.pi/180.)*np.cos(glat*np.pi/180.))/1.e3
yy = adist*np.sin(glon*np.pi/180.)*np.cos(glat*np.pi/180.)/1.e3
zz = adist*np.sin(glat*np.pi/180.)/1.e3
Rgc = np.sqrt(xx**2. + yy**2. + zz**2.)

#open large data set to become greyscale background

b = fits.open('dr17_dr3_McMillan_astroNN_rev1.fits')
hdu = b[1]
agid = b[1].data['GAIAEDR3_SOURCE_ID']
ateff = b[1].data['TEFF']
alogg = b[1].data['LOGG']
asn = b[1].data['SNR']
afe = b[1].data['Fe_H']
amg = b[1].data['Mg_Fe']
ac = b[1].data['C_Fe']
an = b[1].data['N_Fe']
aal = b[1].data['Al_Fe']
amn = b[1].data['Mn_Fe']
asflag = b[1].data['STARFLAG']
adist = b[1].data['weighted_dist']
aedist = b[1].data['weighted_dist_error']
aLz = b[1].data['Lz']
aEnergy = b[1].data['energy']


#find globular clusters that need to be subtracted

gc = fits.open('GC_members_VAC-v1_1.fits')
gcid = gc[1].data['GAIAEDR3_SOURCE_ID']

gcid_mask = np.isin(agid,gcid,invert=True)

#set mask for greyscale background

mask_al = ( (asn >= 50) & (ateff > 3500) & (ateff < 5000) & \
        (alogg < 3.6) & (alogg > -1) & ((aedist/adist)<0.20) & \
        (gcid_mask==True) & (afe > -10) & (amg > -10) & \
        (amn > -10) & (aEnergy < 0.0) & (aLz < 1.e4) & \
        (aLz > -1.e4) & (asflag == 0) & (an < 10) & (an > -10) & \
        (ac > -10) & (ac < 10) & (aal > -10) & (aal < 10))

#Graphs

#work out which stars are from which dwarf galaxy by creating masks for each galaxy
Carina = (field=='CARINA')
Draco = (field=='DRACO')
Fornax = (field=='FORNAX')
Sculptor = (field=='SCULPTOR')
Sextans = (field=='SEXTANS')
Ursa_Minor = (field=='URMINOR')

#Nitrogen Against Iron

fig_size[0] = 64

plt.hist2d(afe[mask_al],an[mask_al],norm=mpl.colors.LogNorm(), bins=(200, 200),range=[[-3, 1],[-1, 2]], cmap='bone')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=tcks)
plt.scatter(feh[Carina], nfe[Carina], c='Aqua', alpha=0.9, s=75, label = 'Carina')
plt.scatter(feh[Draco], nfe[Draco], c='Blue', alpha=0.9, s=75, label = 'Draco')
plt.scatter(feh[Fornax], nfe[Fornax], c='Red', alpha=0.9, s=75, label = 'Fornax')
plt.scatter(feh[Sculptor], nfe[Sculptor], c='Green', alpha=0.9, s=75, label = 'Sculptor')
plt.scatter(feh[Sextans], nfe[Sextans], c='Fuchsia', alpha=0.9, s=75, label = 'Sextans')
plt.scatter(feh[Ursa_Minor], nfe[Ursa_Minor], c='Orange', alpha=0.9, s=75, label = 'Ursa Minor')
plt.ylabel('[N/Fe]', size=labs)
plt.yticks(np.arange(-1, 2, step=1), fontsize=tcks)
plt.xlabel('[Fe/H]', size=labs)
plt.xticks(np.arange(-3, 1, step=1), fontsize=tcks)
plt.gca().set_box_aspect(1)
plt.legend(fontsize = 20, loc='lower right')
plt.savefig('Nitrogen against Iron', bbox_inches='tight')
plt.show()

#Nitrogen against Carbon

fig_size[0] = 64

plt.hist2d(ac[mask_al],an[mask_al],norm=mpl.colors.LogNorm(), bins=(200, 200),range=[[-2, 1.25],[-0.75, 1.5]], cmap='bone')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=tcks)
plt.scatter(cfe[Carina], nfe[Carina], c='Aqua', alpha=0.9, s=75, label = 'Carina')
plt.scatter(cfe[Draco], nfe[Draco], c='Blue', alpha=0.9, s=75, label = 'Draco')
plt.scatter(cfe[Fornax], nfe[Fornax], c='Red', alpha=0.9, s=75, label = 'Fornax')
plt.scatter(cfe[Sculptor], nfe[Sculptor], c='Green', alpha=0.9, s=75, label = 'Sculptor')
plt.scatter(cfe[Sextans], nfe[Sextans], c='Fuchsia', alpha=0.9, s=75, label = 'Sextans')
plt.scatter(cfe[Ursa_Minor], nfe[Ursa_Minor], c='Orange', alpha=0.9, s=75, label = 'Ursa Minor')
plt.ylabel('[N/Fe]', size=labs)
plt.yticks(np.arange(-0.75, 1.5, step=1), fontsize=tcks)
plt.xlabel('[C/Fe]', size=labs)
plt.xticks(np.arange(-2, 1.25, step=1), fontsize=tcks)
plt.gca().set_box_aspect(1)
plt.legend(fontsize = 20, loc='lower right')
plt.savefig('Nitrogen against Carbon', bbox_inches='tight')
plt.show()

#Aluminium against Iron

fig_size[0] = 64

plt.hist2d(afe[mask_al],aal[mask_al],norm=mpl.colors.LogNorm(), bins=(200, 200),range=[[-2.5, 1],[-1.5, 1]], cmap='bone')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=tcks)
plt.scatter(feh[Carina], alfe[Carina], c='Aqua', alpha=0.9, s=75, label = 'Carina')
plt.scatter(feh[Draco], alfe[Draco], c='Blue', alpha=0.9, s=75, label = 'Draco')
plt.scatter(feh[Fornax], alfe[Fornax], c='Red', alpha=0.9, s=75, label = 'Fornax')
plt.scatter(feh[Sculptor], alfe[Sculptor], c='Green', alpha=0.9, s=75, label = 'Sculptor')
plt.scatter(feh[Sextans], alfe[Sextans], c='Fuchsia', alpha=0.9, s=75, label = 'Sextans')
plt.scatter(feh[Ursa_Minor], alfe[Ursa_Minor], c='Orange', alpha=0.9, s=75, label = 'Ursa Minor')
plt.ylabel('[Al/Fe]', size=labs)
plt.yticks(np.arange(-1.5, 1, step=1), fontsize=tcks)
plt.xlabel('[Fe/H]', size=labs)
plt.xticks(np.arange(-2.5, 1, step=1), fontsize=tcks)
plt.gca().set_box_aspect(1)
plt.legend(fontsize = 20, loc='lower right')
plt.savefig('Aluminium against Iron', bbox_inches='tight')
plt.show()

#Magnesium against Aluminium

fig_size[0] = 64

plt.hist2d(amg[mask_al],aal[mask_al],norm=mpl.colors.LogNorm(), bins=(200, 200),range=[[-0.75, 0.75],[-1.5, 1.5]], cmap='bone')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=tcks)
plt.scatter(mgfe[Carina], alfe[Carina], c='Aqua', alpha=0.9, s=75, label = 'Carina')
plt.scatter(mgfe[Draco], alfe[Draco], c='Blue', alpha=0.9, s=75, label = 'Draco')
plt.scatter(mgfe[Fornax], alfe[Fornax], c='Red', alpha=0.9, s=75, label = 'Fornax')
plt.scatter(mgfe[Sculptor], alfe[Sculptor], c='Green', alpha=0.9, s=75, label = 'Sculptor')
plt.scatter(mgfe[Sextans], alfe[Sextans], c='Fuchsia', alpha=0.9, s=75, label = 'Sextans')
plt.scatter(mgfe[Ursa_Minor], alfe[Ursa_Minor], c='Orange', alpha=0.9, s=75, label = 'Ursa Minor')
plt.xlabel('[Mg/Fe]', size=labs)
plt.xticks(np.arange(-0.75, 0.75, step=1), fontsize=tcks)
plt.ylabel('[Al/Fe]', size=labs)
plt.yticks(np.arange(-1.5, 1.5, step=1), fontsize=tcks)
plt.gca().set_box_aspect(1)
plt.legend(fontsize = 20, loc='lower right')
plt.savefig('Magnesium against Aluminium', bbox_inches='tight')
plt.show()