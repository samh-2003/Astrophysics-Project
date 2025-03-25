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
import matplotlib.gridspec as gridspec
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

fig = plt.figure(figsize=(30, 30))
gs = gridspec.GridSpec(2, 2, figure=fig, wspace=0.4, hspace = 0.05)

#Nitrogen Against Iron

ax1 = fig.add_subplot(gs[0, 0])
h1 = ax1.hist2d(afe[mask_al],an[mask_al],norm=mpl.colors.LogNorm(), bins=(200, 200),range=[[-3, 1],[-1, 2]], cmap='bone')
cb1 = fig.colorbar(h1[3], ax=ax1, fraction=0.046, pad=0.04)
cb1.set_label('Counts', fontsize=labs)
cb1.ax.tick_params(labelsize=tcks)
ax1.scatter(feh[Carina], nfe[Carina], c='Aqua', alpha=0.9, s=75, label = 'Carina')
ax1.scatter(feh[Draco], nfe[Draco], c='Blue', alpha=0.9, s=75, label = 'Draco')
ax1.scatter(feh[Fornax], nfe[Fornax], c='Red', alpha=0.9, s=75, label = 'Fornax')
ax1.scatter(feh[Sculptor], nfe[Sculptor], c='Green', alpha=0.9, s=75, label = 'Sculptor')
ax1.scatter(feh[Sextans], nfe[Sextans], c='Fuchsia', alpha=0.9, s=75, label = 'Sextans')
ax1.scatter(feh[Ursa_Minor], nfe[Ursa_Minor], c='Orange', alpha=0.9, s=75, label = 'Ursa Minor')
ax1.set_ylabel('[N/Fe]', size=labs)
ax1.set_yticks(np.arange(-1, 2, step=0.5))
ax1.set_yticklabels(np.arange(-1, 2, step=0.5), fontsize=tcks)
ax1.set_xlabel('[Fe/H]', size=labs)
ax1.set_xticks(np.arange(-3, 1, step=0.5))
ax1.set_xticklabels(np.arange(-3, 1, step=0.5), fontsize=tcks, rotation = 45)
ax1.set_box_aspect(1)
ax1.legend(fontsize = 20, loc='lower right')


#Nitrogen against Carbon

ax2 = fig.add_subplot(gs[0, 1])
h2 = ax2.hist2d(ac[mask_al],an[mask_al],norm=mpl.colors.LogNorm(), bins=(200, 200),range=[[-2, 1.25],[-0.75, 1.5]], cmap='bone')
cb2 = fig.colorbar(h2[3], ax=ax2, fraction=0.046, pad=0.04)
cb2.set_label('Counts', fontsize=labs)
cb2.ax.tick_params(labelsize=tcks)
ax2.scatter(cfe[Carina], nfe[Carina], c='Aqua', alpha=0.9, s=75, label = 'Carina')
ax2.scatter(cfe[Draco], nfe[Draco], c='Blue', alpha=0.9, s=75, label = 'Draco')
ax2.scatter(cfe[Fornax], nfe[Fornax], c='Red', alpha=0.9, s=75, label = 'Fornax')
ax2.scatter(cfe[Sculptor], nfe[Sculptor], c='Green', alpha=0.9, s=75, label = 'Sculptor')
ax2.scatter(cfe[Sextans], nfe[Sextans], c='Fuchsia', alpha=0.9, s=75, label = 'Sextans')
ax2.scatter(cfe[Ursa_Minor], nfe[Ursa_Minor], c='Orange', alpha=0.9, s=75, label = 'Ursa Minor')
ax2.set_ylabel('[N/Fe]', size=labs)
ax2.set_yticks(np.arange(-0.75, 1.5, step=0.5))
ax2.set_yticklabels(np.arange(-0.75, 1.5, step=0.5), fontsize=tcks)
ax2.set_xlabel('[C/Fe]', size=labs)
ax2.set_xticks(np.arange(-2, 1.25, step=0.5))
ax2.set_xticklabels(np.arange(-2, 1.25, step=0.5), fontsize=tcks, rotation = 45)
ax2.set_box_aspect(1)
ax2.legend(fontsize = 20, loc='lower right')

#Aluminium against Iron

ax3 = fig.add_subplot(gs[1, 0])
h3 = ax3.hist2d(afe[mask_al],aal[mask_al],norm=mpl.colors.LogNorm(), bins=(200, 200),range=[[-2.5, 1],[-1.5, 1]], cmap='bone')
cb3 = fig.colorbar(h3[3], ax=ax3, fraction=0.046, pad=0.04)
cb3.set_label('Counts', fontsize=labs)
cb3.ax.tick_params(labelsize=tcks)
ax3.scatter(feh[Carina], alfe[Carina], c='Aqua', alpha=0.9, s=75, label = 'Carina')
ax3.scatter(feh[Draco], alfe[Draco], c='Blue', alpha=0.9, s=75, label = 'Draco')
ax3.scatter(feh[Fornax], alfe[Fornax], c='Red', alpha=0.9, s=75, label = 'Fornax')
ax3.scatter(feh[Sculptor], alfe[Sculptor], c='Green', alpha=0.9, s=75, label = 'Sculptor')
ax3.scatter(feh[Sextans], alfe[Sextans], c='Fuchsia', alpha=0.9, s=75, label = 'Sextans')
ax3.scatter(feh[Ursa_Minor], alfe[Ursa_Minor], c='Orange', alpha=0.9, s=75, label = 'Ursa Minor')
ax3.set_ylabel('[Al/Fe]', size=labs)
ax3.set_yticks(np.arange(-1.5, 1, step=0.5))
ax3.set_yticklabels(np.arange(-1.5, 1, step=0.5), fontsize=tcks)
ax3.set_xlabel('[Fe/H]', size=labs)
ax3.set_xticks(np.arange(-2.5, 1, step=0.5))
ax3.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks, rotation = 45)
ax3.set_box_aspect(1)
ax3.legend(fontsize = 20, loc='lower right')

#Magnesium against Aluminium

ax4 = fig.add_subplot(gs[1, 1])
h4 = ax4.hist2d(amg[mask_al],aal[mask_al],norm=mpl.colors.LogNorm(), bins=(200, 200),range=[[-0.75, 0.75],[-1.5, 1.5]], cmap='bone')
cb4 = fig.colorbar(h4[3], ax=ax4, fraction=0.046, pad=0.04)
cb4.set_label('Counts', fontsize=labs)
cb4.ax.tick_params(labelsize=tcks)
ax4.scatter(mgfe[Carina], alfe[Carina], c='Aqua', alpha=0.9, s=75, label = 'Carina')
ax4.scatter(mgfe[Draco], alfe[Draco], c='Blue', alpha=0.9, s=75, label = 'Draco')
ax4.scatter(mgfe[Fornax], alfe[Fornax], c='Red', alpha=0.9, s=75, label = 'Fornax')
ax4.scatter(mgfe[Sculptor], alfe[Sculptor], c='Green', alpha=0.9, s=75, label = 'Sculptor')
ax4.scatter(mgfe[Sextans], alfe[Sextans], c='Fuchsia', alpha=0.9, s=75, label = 'Sextans')
ax4.scatter(mgfe[Ursa_Minor], alfe[Ursa_Minor], c='Orange', alpha=0.9, s=75, label = 'Ursa Minor')
ax4.set_xlabel('[Mg/Fe]', size=labs)
ax4.set_xticks(np.arange(-0.75, 0.75, step=0.25))
ax4.set_xticklabels(np.arange(-0.75, 0.75, step=0.25), fontsize=tcks, rotation = 45)
ax4.set_ylabel('[Al/Fe]', size=labs)
ax4.set_yticks(np.arange(-1.5, 1.5, step=0.5))
ax4.set_yticklabels(np.arange(-1.5, 1.5, step=0.5), fontsize=tcks)
ax4.set_box_aspect(1)
ax4.legend(fontsize = 20, loc='lower right')


plt.show()