# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:06:41 2025

@author: samhe
"""

#Abundance graphs
from astropy.table import Table
from astropy.io import fits 
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import pylab as plt 
import numpy as np
import math
from scipy import stats
import pandas as pd 
from pandas import DataFrame

# Get current size
fig_size = plt.rcParams["figure.figsize"]

# Set figure width to 16 and height to 12
fig_size[0] = 20
fig_size[1] = 20

#Figure settings (font sizes of ticks and labels)
labs=40
tcks=35

#Open G2 Files
G2 = fits.open('DwarfG2s.fits')

ids_g2 = G2[1].data['APOGEE_ID_1']
afield_g2 = G2[1].data['FIELD']
mgfe_g2 = G2[1].data['Mg_Fe']
nfe_g2 = G2[1].data['N_Fe']
nife_g2 = G2[1].data['Ni_Fe']
feh_g2 = G2[1].data['Fe_H']
cfe_g2 = G2[1].data['C_Fe']
cafe_g2 = G2[1].data['CA_Fe']
sife_g2 = G2[1].data['Si_Fe']
alfe_g2 = G2[1].data['Al_Fe']
mnfe_g2 = G2[1].data['Mn_Fe']


#Open Dwarf galaxy files
DG = fits.open('DwarfGalaxies.fits')

ids_al = DG[1].data['APOGEE_ID_1']
afield_al = DG[1].data['FIELD']
mgfe_al = DG[1].data['Mg_Fe']
nfe_al = DG[1].data['N_Fe']
nife_al = DG[1].data['Ni_Fe']
feh_al = DG[1].data['Fe_H']
cfe_al = DG[1].data['C_Fe']
cafe_al = DG[1].data['CA_Fe']
sife_al = DG[1].data['Si_Fe']
alfe_al = DG[1].data['Al_Fe']
mnfe_al = DG[1].data['Mn_Fe']

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
ani = b[1].data['Ni_Fe']
asi = b[1].data['Si_Fe']
aca = b[1].data['CA_FE']
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

gs = gridspec.GridSpec(2,5)
gs.update()
#wspace=0.2, hspace=0.1

#Mg against Fe
plt.subplot(gs[0,0])
plt.rc('text', usetex=True)
plt.suptitle(r'$\underline{[Fe/H] against Metallicities}$',y = 0.94, fontsize = 50, usetex = True)
plt.hist2d(afe[mask_al],amg[mask_al],norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
plt.scatter(feh_al, mgfe_al, c='green', alpha=1.0, s=50, label = 'Dwarf Galaxies')
plt.scatter(feh_g2, mgfe_g2, c='red', alpha=0.9, s=100, label = 'G2s')
plt.xlabel('[Fe/H]', size=labs)
plt.xticks(np.arange(-2.5, 1, step=0.5), fontsize=tcks)
plt.ylabel('[Mg/Fe]', size=labs)
plt.yticks(np.arange(-0.75, 0.75, step=0.25), fontsize=tcks)
plt.gca().set_box_aspect(1)
plt.legend(fontsize = 20, loc='lower right')

#Al against Fe
plt.subplot(gs[0,1])
plt.hist2d(afe[mask_al],aal[mask_al],norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
plt.scatter(feh_al, alfe_al, c='green', alpha=1.0, s=50, label = 'Dwarf Galaxies')
plt.scatter(feh_g2, alfe_g2, c='red', alpha=0.9, s=100, label = 'G2s')
plt.xlabel('[Fe/H]', size=labs)
plt.xticks(np.arange(-2.5, 1, step=0.5), fontsize=tcks)
plt.ylabel('[Al/Fe]', size=labs)
plt.yticks(np.arange(-0.75, 0.75, step=0.25), fontsize=tcks)
plt.gca().set_box_aspect(1)
plt.legend(fontsize = 20, loc='lower right')


#Si against Fe
plt.subplot(gs[0,2])
plt.hist2d(afe[mask_al],asi[mask_al],norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
plt.scatter(feh_al, sife_al, c='green', alpha=1.0, s=50, label = 'Dwarf Galaxies')
plt.scatter(feh_g2, sife_g2, c='red', alpha=0.9, s=100, label = 'G2s')
plt.xlabel('[Fe/H]', size=labs)
plt.xticks(np.arange(-2.5, 1, step=0.5), fontsize=tcks)
plt.ylabel('[Si/Fe]', size=labs)
plt.yticks(np.arange(-0.75, 0.75, step=0.25), fontsize=tcks)
plt.gca().set_box_aspect(1)
plt.legend(fontsize = 20, loc='lower right')

#Ni against Fe
plt.subplot(gs[0,3])
plt.rc('text', usetex=True)
plt.suptitle(r'$\underline{[Fe/H] against Metallicities}$',y = 0.94, fontsize = 50, usetex = True)
#plt.hist2d(afe[mask_al],ani[mask_al],norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
plt.scatter(feh_al, nife_al, c='green', alpha=1.0, s=50, label = 'Dwarf Galaxies')
plt.scatter(feh_g2, nife_g2, c='red', alpha=0.9, s=100, label = 'G2s')
plt.xlabel('[Fe/H]', size=labs)
plt.xticks(np.arange(-2.5, 1, step=0.5), fontsize=tcks)
plt.ylabel('[Ni/Fe]', size=labs)
plt.yticks(np.arange(-0.75, 0.75, step=0.25), fontsize=tcks)
plt.gca().set_box_aspect(1)
plt.legend(fontsize = 20, loc='lower right')

#C against Fe
plt.subplot(gs[1,0])
#plt.hist2d(afe[mask_al],ac[mask_al],norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
plt.scatter(feh_al, cfe_al, c='green', alpha=1.0, s=50, label = 'Dwarf Galaxies')
plt.scatter(feh_g2, cfe_g2, c='red', alpha=0.9, s=100, label = 'G2s')
plt.xlabel('[Fe/H]', size=labs)
plt.xticks(np.arange(-2.5, 1, step=0.5), fontsize=tcks)
plt.ylabel('[C/Fe]', size=labs)
plt.yticks(np.arange(-0.75, 0.75, step=0.25), fontsize=tcks)
plt.gca().set_box_aspect(1)
plt.legend(fontsize = 20, loc='lower right')

#N against Fe
plt.subplot(gs[1,1])
#plt.hist2d(afe[mask_al],an[mask_al],norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
plt.scatter(feh_al, nfe_al, c='green', alpha=1.0, s=50, label = 'Dwarf Galaxies')
plt.scatter(feh_g2, nfe_g2, c='red', alpha=0.9, s=100, label = 'G2s')
plt.xlabel('[Fe/H]', size=labs)
plt.xticks(np.arange(-2.5, 1, step=0.5), fontsize=tcks)
plt.ylabel('[N/Fe]', size=labs)
plt.yticks(np.arange(-0.75, 0.75, step=0.25), fontsize=tcks)
plt.gca().set_box_aspect(1)
plt.legend(fontsize = 20, loc='lower right')

#Mn against Fe
plt.subplot(gs[1,2])
#plt.hist2d(afe[mask_al],amn[mask_al],norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
plt.scatter(feh_al, mnfe_al, c='green', alpha=1.0, s=50, label = 'Dwarf Galaxies')
plt.scatter(feh_g2, mnfe_g2, c='red', alpha=0.9, s=100, label = 'G2s')
plt.xlabel('[Fe/H]', size=labs)
plt.xticks(np.arange(-2.5, 1, step=0.5), fontsize=tcks)
plt.ylabel('[Mn/Fe]', size=labs)
plt.yticks(np.arange(-0.75, 0.75, step=0.25), fontsize=tcks)
plt.gca().set_box_aspect(1)
plt.legend(fontsize = 20, loc='lower right')

#Ca against Fe
plt.subplot(gs[1,3])
#plt.hist2d(afe[mask_al],aca[mask_al],norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
plt.scatter(feh_al, cafe_al, c='green', alpha=1.0, s=50, label = 'Dwarf Galaxies')
plt.scatter(feh_g2, cafe_g2, c='red', alpha=0.9, s=100, label = 'G2s')
plt.xlabel('[Fe/H]', size=labs)
plt.xticks(np.arange(-2.5, 1, step=0.5), fontsize=tcks)
plt.ylabel('[Ca/Fe]', size=labs)
plt.yticks(np.arange(-0.75, 0.75, step=0.25), fontsize=tcks)
plt.gca().set_box_aspect(1)
plt.legend(fontsize = 20, loc='lower right')

ax = plt.subplot(gs[:,4])
cb = plt.colorbar(ax)
cb.ax.tick_params(labelsize=tcks)