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
import numpy as np
import math
from scipy import stats
import pandas as pd 
from pandas import DataFrame

# Get current size
fig_size = plt.rcParams["figure.figsize"]

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
fig = plt.figure(figsize=(30, 30))
fig2 = plt.figure(figsize=(30, 30))
gs = gridspec.GridSpec(2, 2, figure=fig, wspace=0.4, hspace=0.1)
gs2 = gridspec.GridSpec(2, 2, figure=fig2, wspace=0.4, hspace=0.1)


# Mg against Fe
ax1 = fig.add_subplot(gs[0, 0])
h1 = ax1.hist2d(afe[mask_al], amg[mask_al], norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb1 = fig.colorbar(h1[3], ax=ax1, fraction=0.046, pad=0.04)
cb1.set_label('Counts', fontsize=labs)
cb1.ax.tick_params(labelsize=tcks)
ax1.scatter(feh_al, mgfe_al, c='green', alpha=1.0, s=50, label='Dwarf Galaxies')
ax1.scatter(feh_g2, mgfe_g2, c='red', alpha=0.9, s=100, label='G2s')
ax1.set_xlabel('[Fe/H]', size=labs)
ax1.set_xticks(np.arange(-2.5, 1, step=0.5))
ax1.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks)
ax1.set_ylabel('[Mg/Fe]', size=labs)
ax1.set_yticks(np.arange(-0.75, 0.75, step=0.25))
ax1.set_yticklabels(np.arange(-0.75, 0.75, step=0.25), fontsize=tcks)
ax1.set_box_aspect(1)
ax1.legend(fontsize=20, loc='lower right')

# Al against Fe
ax2 = fig.add_subplot(gs[0, 1])
h2 = ax2.hist2d(afe[mask_al], aal[mask_al], norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb2 = fig.colorbar(h2[3], ax=ax2, fraction=0.046, pad=0.04)
cb2.set_label('Counts', fontsize=labs)
cb2.ax.tick_params(labelsize=tcks)
ax2.scatter(feh_al, alfe_al, c='green', alpha=1.0, s=50, label='Dwarf Galaxies')
ax2.scatter(feh_g2, alfe_g2, c='red', alpha=0.9, s=100, label='G2s')
ax2.set_xlabel('[Fe/H]', size=labs)
ax2.set_xticks(np.arange(-2.5, 1, step=0.5))
ax2.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks)
ax2.set_ylabel('[Al/Fe]', size=labs)
ax2.set_yticks(np.arange(-1, 1, step=0.25))
ax2.set_yticklabels(np.arange(-1, 1, step=0.25), fontsize=tcks)
ax2.set_box_aspect(1)
ax2.legend(fontsize=20, loc='lower right')

# Si against Fe
ax3 = fig.add_subplot(gs[1, 0])
h3 = ax3.hist2d(afe[mask_al], asi[mask_al], norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb3 = fig.colorbar(h3[3], ax=ax3, fraction=0.046, pad=0.04)
cb3.set_label('Counts', fontsize=labs)
cb3.ax.tick_params(labelsize=tcks)
ax3.scatter(feh_al, sife_al, c='green', alpha=1.0, s=50, label='Dwarf Galaxies')
ax3.scatter(feh_g2, sife_g2, c='red', alpha=0.9, s=100, label='G2s')
ax3.set_xlabel('[Fe/H]', size=labs)
ax3.set_xticks(np.arange(-2.5, 1, step=0.5))
ax3.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks)
ax3.set_ylabel('[Si/Fe]', size=labs)
ax3.set_yticks(np.arange(-0.75, 0.75, step=0.25))
ax3.set_yticklabels(np.arange(-0.75, 0.75, step=0.25), fontsize=tcks)
ax3.set_box_aspect(1)
ax3.legend(fontsize=20, loc='lower right')

# Ni against Fe

#replace nan values from ani with 0
ani_clean = np.nan_to_num(ani[mask_al])
ax4 = fig.add_subplot(gs[1, 1])
h4 = ax4.hist2d(afe[mask_al], ani_clean, norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb4 = fig.colorbar(h4[3], ax=ax4, fraction=0.046, pad=0.04)
cb4.set_label('Counts', fontsize=labs)
cb4.ax.tick_params(labelsize=tcks)
ax4.scatter(feh_al, nife_al, c='green', alpha=1.0, s=50, label='Dwarf Galaxies')
ax4.scatter(feh_g2, nife_g2, c='red', alpha=0.9, s=100, label='G2s')
ax4.set_xlabel('[Fe/H]', size=labs)
ax4.set_xticks(np.arange(-2.5, 1, step=0.5))
ax4.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks)
ax4.set_ylabel('[Ni/Fe]', size=labs)
ax4.set_yticks(np.arange(-0.75, 0.75, step=0.25))
ax4.set_yticklabels(np.arange(-0.75, 0.75, step=0.25), fontsize=tcks)
ax4.set_box_aspect(1)
ax4.legend(fontsize=20, loc='lower right')


# C against Fe
ax5 = fig2.add_subplot(gs2[0, 0])
h5 = ax5.hist2d(afe[mask_al], ac[mask_al], norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb5 = fig2.colorbar(h5[3], ax=ax5, fraction=0.046, pad=0.04)
cb5.set_label('Counts', fontsize=labs)
cb5.ax.tick_params(labelsize=tcks)
ax5.scatter(feh_al, cfe_al, c='green', alpha=1.0, s=50, label='Dwarf Galaxies')
ax5.scatter(feh_g2, cfe_g2, c='red', alpha=0.9, s=100, label='G2s')
ax5.set_xlabel('[Fe/H]', size=labs)
ax5.set_xticks(np.arange(-2.5, 1, step=0.5))
ax5.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks)
ax5.set_ylabel('[C/Fe]', size=labs)
ax5.set_yticks(np.arange(-1.5, 1, step=0.5))
ax5.set_yticklabels(np.arange(-1.5, 1, step=0.5), fontsize=tcks)
ax5.set_box_aspect(1)
ax5.legend(fontsize=20, loc='lower right')

# N against Fe
ax6 = fig2.add_subplot(gs2[0, 1])
h6 = ax6.hist2d(afe[mask_al], an[mask_al], norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb6 = fig2.colorbar(h6[3], ax=ax6, fraction=0.046, pad=0.04)
cb6.set_label('Counts', fontsize=labs)
cb6.ax.tick_params(labelsize=tcks)
ax6.scatter(feh_al, nfe_al, c='green', alpha=1.0, s=50, label='Dwarf Galaxies')
ax6.scatter(feh_g2, nfe_g2, c='red', alpha=0.9, s=100, label='G2s')
ax6.set_xlabel('[Fe/H]', size=labs)
ax6.set_xticks(np.arange(-2.5, 1, step=0.5))
ax6.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks)
ax6.set_ylabel('[N/Fe]', size=labs)
ax6.set_yticks(np.arange(-0.75, 2, step=0.5))
ax6.set_yticklabels(np.arange(-0.75, 2, step=0.5), fontsize=tcks)
ax6.set_box_aspect(1)
ax6.legend(fontsize=20, loc='lower right')

# Mn against Fe
ax7 = fig2.add_subplot(gs2[1, 0])
h7 = ax7.hist2d(afe[mask_al], amn[mask_al], norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb7 = fig2.colorbar(h7[3], ax=ax7, fraction=0.046, pad=0.04)
cb7.set_label('Counts', fontsize=labs)
cb7.ax.tick_params(labelsize=tcks)
ax7.scatter(feh_al, mnfe_al, c='green', alpha=1.0, s=50, label='Dwarf Galaxies')
ax7.scatter(feh_g2, mnfe_g2, c='red', alpha=0.9, s=100, label='G2s')
ax7.set_xlabel('[Fe/H]', size=labs)
ax7.set_xticks(np.arange(-2.5, 1, step=0.5))
ax7.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks)
ax7.set_ylabel('[Mn/Fe]', size=labs)
ax7.set_yticks(np.arange(-1, 1.5, step=0.5))
ax7.set_yticklabels(np.arange(-1, 1.5, step=0.5), fontsize=tcks)
ax7.set_box_aspect(1)
ax7.legend(fontsize=20, loc='lower right')

# Ca against Fe
#replace nan values from aca with 0
aca_clean = np.nan_to_num(aca[mask_al])

ax8 = fig2.add_subplot(gs2[1, 1])
h8 = ax8.hist2d(afe[mask_al], aca_clean, norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb8 = fig2.colorbar(h8[3], ax=ax8, fraction=0.046, pad=0.04)
cb8.set_label('Counts', fontsize=labs)
cb8.ax.tick_params(labelsize=tcks)
ax8.scatter(feh_al, cafe_al, c='green', alpha=1.0, s=50, label='Dwarf Galaxies')
ax8.scatter(feh_g2, cafe_g2, c='red', alpha=0.9, s=100, label='G2s')
ax8.set_xlabel('[Fe/H]', size=labs)
ax8.set_xticks(np.arange(-2.5, 1, step=0.5))
ax8.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks)
ax8.set_ylabel('[Ca/Fe]', size=labs)
ax8.set_yticks(np.arange(-0.75, 1, step=0.25))
ax8.set_yticklabels(np.arange(-0.75, 1, step=0.25), fontsize=tcks)
ax8.set_box_aspect(1)
ax8.legend(fontsize=20, loc='lower right')

plt.show()