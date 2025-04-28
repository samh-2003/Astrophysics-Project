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
from matplotlib.lines import Line2D
from matplotlib.container import ErrorbarContainer

# Get current size
fig_size = plt.rcParams["figure.figsize"]

#Figure settings (font sizes of ticks and labels)
labs=30
tcks=25
legend=17.5

#Open G2 File
G2 = fits.open('NewDwarfG2s.fits')

hdu = G2[1]
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
ofe_g2 = G2[1].data['O_Fe']
errmgfe_g2 = G2[1].data['Mg_Fe_Err']
errnfe_g2 = G2[1].data['N_Fe_Err']
errnife_g2 = G2[1].data['Ni_Fe_Err']
errfeh_g2 = G2[1].data['Fe_H_Err']
errcfe_g2 = G2[1].data['C_Fe_Err']
errcafe_g2 = G2[1].data['CA_Fe_Err']
errsife_g2 = G2[1].data['Si_Fe_Err']
erralfe_g2 = G2[1].data['Al_Fe_Err']
errmnfe_g2 = G2[1].data['Mn_Fe_Err']
errofe_g2 = G2[1].data['O_Fe_Err']
print(afield_g2)

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
ofe_al = DG[1].data['O_Fe']

#open large data set to become greyscale background

b = fits.open('dr17_dr3_McMillan_astroNN_rev1.fits')
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
ao = b[1].data['O_Fe']
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

# Function to calculate binned mean of x and y errors
def bin_mean_errors(x, y, x_err, y_err, x_bins):
    bin_centers = []
    bin_x_error_means = []  # Mean of x errors
    bin_y_error_means = []  # Mean of y errors
    for i in range(len(x_bins) - 1):
        mask = (x >= x_bins[i]) & (x < x_bins[i + 1])
        if np.sum(mask) > 0:
            bin_center = (x_bins[i] + x_bins[i + 1]) / 2
            bin_centers.append(bin_center)
            bin_x_error_means.append(np.nanmean(x_err[mask]))  # Mean of x errors
            bin_y_error_means.append(np.nanmean(y_err[mask]))  # Mean of y errors
    bin_x_error_means = np.array(bin_x_error_means, dtype=float).flatten()
    bin_y_error_means = np.array(bin_y_error_means, dtype=float).flatten()
    return np.array(bin_centers), np.array(bin_x_error_means), np.array(bin_y_error_means)

# Define bins
x_bins = np.arange(-2.5, -0.25, 0.25)

#define masks for each dwarf galaxy, to identify which G2 object is in what dwarf galaxy
Carina = (afield_g2=='CARINA')
Draco = (afield_g2=='DRACO')
Fornax = (afield_g2=='FORNAX')
Sculptor = (afield_g2=='SCULPTOR')
Ursa_Minor = (afield_g2=='URMINOR')
#sextans has no G2 stars, so no need to code for it


#Graphs
fig = plt.figure(figsize=(30, 30))
gs = gridspec.GridSpec(3, 3, figure=fig, wspace=0.375, hspace=0)

plt.rc('text', usetex=True)

# Mg against Fe
ax1 = fig.add_subplot(gs[0, 0])
h1 = ax1.hist2d(afe[mask_al], amg[mask_al], norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb1 = fig.colorbar(h1[3], ax=ax1, fraction=0.046, pad=0.04)
cb1.set_label('Counts', fontsize=labs)
cb1.ax.tick_params(labelsize=tcks)
bin_centers, bin_x_error_means, bin_y_error_means = bin_mean_errors(feh_g2, mgfe_g2, errfeh_g2, errmgfe_g2, x_bins)
ax1.errorbar(bin_centers, [-0.65] * len(bin_centers), xerr=bin_x_error_means, yerr=bin_y_error_means, 
             ls = 'none', color='blue', label='Mean Errors')
ax1.scatter(feh_al, mgfe_al, c='green', alpha=1.0, s=15, label='Dwarf Galaxies')
ax1.scatter(feh_g2[Carina], mgfe_g2[Carina], c='red', alpha=0.9, s=50, label='Carina G2s')
ax1.scatter(feh_g2[Draco], mgfe_g2[Draco], c='m', alpha=0.9, s=50, label='Draco G2s')
ax1.scatter(feh_g2[Fornax], mgfe_g2[Fornax], c='orange', alpha=0.9, s=50, label='Fornax G2s')
ax1.scatter(feh_g2[Sculptor], mgfe_g2[Sculptor], c='cyan', alpha=0.9, s=50, label='Sculptor G2s')
ax1.scatter(feh_g2[Ursa_Minor], mgfe_g2[Ursa_Minor], c='dodgerblue', alpha=0.9, s=50, label='Ursa Minor G2s')
ax1.set_xlabel('[Fe/H]', size=labs)
ax1.set_xticks(np.arange(-2.5, 1, step=0.5))
ax1.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks, rotation = 45)
ax1.set_ylabel('[Mg/Fe]', size=labs)
ax1.set_yticks(np.arange(-0.75, 0.75, step=0.25))
ax1.set_yticklabels(np.arange(-0.75, 0.75, step=0.25), fontsize=tcks)
ax1.set_box_aspect(1)
#code to remove the marker from the errorbar
handles, labels = ax1.get_legend_handles_labels()
handles = [h[0] if isinstance(h, ErrorbarContainer) else h for h in handles]
ax1.legend(handles, labels, fontsize=legend, loc='lower right', frameon=False, labelcolor='linecolor', markerscale=0, labelspacing = 0.1)

# Al against Fe
ax2 = fig.add_subplot(gs[0, 1])
h2 = ax2.hist2d(afe[mask_al], aal[mask_al], norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb2 = fig.colorbar(h2[3], ax=ax2, fraction=0.046, pad=0.04)
cb2.set_label('Counts', fontsize=labs)
cb2.ax.tick_params(labelsize=tcks)
bin_centers, bin_x_error_means, bin_y_error_means = bin_mean_errors(feh_g2, alfe_g2, errfeh_g2, erralfe_g2, x_bins)
ax2.errorbar(bin_centers, [0.75] * len(bin_centers), xerr=bin_x_error_means, yerr=bin_y_error_means, 
             ls = 'none', color='blue', label='Mean Errors')
ax2.scatter(feh_al, alfe_al, c='green', alpha=1.0, s=15, label='Dwarf Galaxies')
ax2.scatter(feh_g2[Carina], alfe_g2[Carina], c='red', alpha=0.9, s=50, label='Carina G2s')
ax2.scatter(feh_g2[Draco], alfe_g2[Draco], c='m', alpha=0.9, s=50, label='Draco G2s')
ax2.scatter(feh_g2[Fornax], alfe_g2[Fornax], c='orange', alpha=0.9, s=50, label='Fornax G2s')
ax2.scatter(feh_g2[Sculptor], alfe_g2[Sculptor], c='cyan', alpha=0.9, s=50, label='Sculptor G2s')
ax2.scatter(feh_g2[Ursa_Minor], alfe_g2[Ursa_Minor], c='dodgerblue', alpha=0.9, s=50, label='Ursa Minor G2s')
ax2.set_xlabel('[Fe/H]', size=labs)
ax2.set_xticks(np.arange(-2.5, 1, step=0.5))
ax2.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks, rotation = 45)
ax2.set_ylabel('[Al/Fe]', size=labs)
ax2.set_yticks(np.arange(-1, 1, step=0.25))
ax2.set_yticklabels(np.arange(-1, 1, step=0.25), fontsize=tcks)
ax2.set_box_aspect(1)
#code to remove the marker from the errorbar
handles, labels = ax2.get_legend_handles_labels()
handles = [h[0] if isinstance(h, ErrorbarContainer) else h for h in handles]
ax2.legend(handles, labels, fontsize=legend, loc='upper right', frameon=False, labelcolor='linecolor', markerscale=0, labelspacing = 0.1)

# Si against Fe
ax3 = fig.add_subplot(gs[0, 2])
h3 = ax3.hist2d(afe[mask_al], asi[mask_al], norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb3 = fig.colorbar(h3[3], ax=ax3, fraction=0.046, pad=0.04)
cb3.set_label('Counts', fontsize=labs)
cb3.ax.tick_params(labelsize=tcks)
bin_centers, bin_x_error_means, bin_y_error_means = bin_mean_errors(feh_g2, sife_g2, errfeh_g2, errsife_g2, x_bins)
ax3.errorbar(bin_centers, [0.6] * len(bin_centers), xerr=bin_x_error_means, yerr=bin_y_error_means, 
             ls = 'none', color='blue', label='Mean Errors')
ax3.scatter(feh_al, sife_al, c='green', alpha=1.0, s=15, label='Dwarf Galaxies')
ax3.scatter(feh_g2[Carina], sife_g2[Carina], c='red', alpha=0.9, s=50, label='Carina G2s')
ax3.scatter(feh_g2[Draco], sife_g2[Draco], c='m', alpha=0.9, s=50, label='Draco G2s')
ax3.scatter(feh_g2[Fornax], sife_g2[Fornax], c='orange', alpha=0.9, s=50, label='Fornax G2s')
ax3.scatter(feh_g2[Sculptor], sife_g2[Sculptor], c='cyan', alpha=0.9, s=50, label='Sculptor G2s')
ax3.scatter(feh_g2[Ursa_Minor], sife_g2[Ursa_Minor], c='dodgerblue', alpha=0.9, s=50, label='Ursa Minor G2s')
ax3.set_xlabel('[Fe/H]', size=labs)
ax3.set_xticks(np.arange(-2.5, 1, step=0.5))
ax3.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks, rotation = 45)
ax3.set_ylabel('[Si/Fe]', size=labs)
ax3.set_yticks(np.arange(-0.75, 0.75, step=0.25))
ax3.set_yticklabels(np.arange(-0.75, 0.75, step=0.25), fontsize=tcks)
ax3.set_box_aspect(1)
#code to remove the marker from the errorbar
handles, labels = ax3.get_legend_handles_labels()
handles = [h[0] if isinstance(h, ErrorbarContainer) else h for h in handles]
ax3.legend(handles, labels, fontsize=legend, loc='upper right', frameon=False, labelcolor='linecolor', markerscale=0, labelspacing = 0.1)

# Ni against Fe

#replace nan values from ani with 0
ani_clean = np.nan_to_num(ani[mask_al])
ax4 = fig.add_subplot(gs[1, 0])
h4 = ax4.hist2d(afe[mask_al], ani_clean, norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb4 = fig.colorbar(h4[3], ax=ax4, fraction=0.046, pad=0.04)
cb4.set_label('Counts', fontsize=labs)
cb4.ax.tick_params(labelsize=tcks)
bin_centers, bin_x_error_means, bin_y_error_means = bin_mean_errors(feh_g2, nife_g2, errfeh_g2, errnife_g2, x_bins)
ax4.errorbar(bin_centers, [-0.75] * len(bin_centers), xerr=bin_x_error_means, yerr=bin_y_error_means, 
             ls = 'none', color='blue', label='Mean Errors')
ax4.scatter(feh_al, nife_al, c='green', alpha=1.0, s=15, label='Dwarf Galaxies')
ax4.scatter(feh_g2[Carina], nife_g2[Carina], c='red', alpha=0.9, s=50, label='Carina G2s')
ax4.scatter(feh_g2[Draco], nife_g2[Draco], c='m', alpha=0.9, s=50, label='Draco G2s')
ax4.scatter(feh_g2[Fornax], nife_g2[Fornax], c='orange', alpha=0.9, s=50, label='Fornax G2s')
ax4.scatter(feh_g2[Sculptor], nife_g2[Sculptor], c='cyan', alpha=0.9, s=50, label='Sculptor G2s')
ax4.scatter(feh_g2[Ursa_Minor], nife_g2[Ursa_Minor], c='dodgerblue', alpha=0.9, s=50, label='Ursa Minor G2s')
ax4.set_xlabel('[Fe/H]', size=labs)
ax4.set_xticks(np.arange(-2.5, 1, step=0.5))
ax4.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks, rotation = 45)
ax4.set_ylabel('[Ni/Fe]', size=labs)
ax4.set_yticks(np.arange(-0.75, 0.75, step=0.25))
ax4.set_yticklabels(np.arange(-0.75, 0.75, step=0.25), fontsize=tcks)
ax4.set_box_aspect(1)
#code to remove the marker from the errorbar
handles, labels = ax4.get_legend_handles_labels()
handles = [h[0] if isinstance(h, ErrorbarContainer) else h for h in handles]
ax4.legend(handles, labels, fontsize=legend, loc='lower right', frameon=False, labelcolor='linecolor', markerscale=0, labelspacing = 0.1)


# C against Fe
ax5 = fig.add_subplot(gs[1, 1])
h5 = ax5.hist2d(afe[mask_al], ac[mask_al], norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb5 = fig.colorbar(h5[3], ax=ax5, fraction=0.046, pad=0.04)
cb5.set_label('Counts', fontsize=labs)
cb5.ax.tick_params(labelsize=tcks)
bin_centers, bin_x_error_means, bin_y_error_means = bin_mean_errors(feh_g2, cfe_g2, errfeh_g2, errcfe_g2, x_bins)
ax5.errorbar(bin_centers, [0.58] * len(bin_centers), xerr=bin_x_error_means, yerr=bin_y_error_means, 
             ls = 'none', color='blue', label='Mean Errors')
ax5.scatter(feh_al, cfe_al, c='green', alpha=1.0, s=15, label='Dwarf Galaxies')
ax5.scatter(feh_g2[Carina], cfe_g2[Carina], c='red', alpha=0.9, s=50, label='Carina G2s')
ax5.scatter(feh_g2[Draco], cfe_g2[Draco], c='m', alpha=0.9, s=50, label='Draco G2s')
ax5.scatter(feh_g2[Fornax], cfe_g2[Fornax], c='orange', alpha=0.9, s=50, label='Fornax G2s')
ax5.scatter(feh_g2[Sculptor], cfe_g2[Sculptor], c='cyan', alpha=0.9, s=50, label='Sculptor G2s')
ax5.scatter(feh_g2[Ursa_Minor], cfe_g2[Ursa_Minor], c='dodgerblue', alpha=0.9, s=50, label='Ursa Minor G2s')
ax5.set_xlabel('[Fe/H]', size=labs)
ax5.set_xticks(np.arange(-2.5, 1, step=0.5))
ax5.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks, rotation = 45)
ax5.set_ylabel('[C/Fe]', size=labs)
ax5.set_yticks(np.arange(-1.5, 1, step=0.5))
ax5.set_yticklabels(np.arange(-1.5, 1, step=0.5), fontsize=tcks)
ax5.set_box_aspect(1)
#code to remove the marker from the errorbar
handles, labels = ax5.get_legend_handles_labels()
handles = [h[0] if isinstance(h, ErrorbarContainer) else h for h in handles]
ax5.legend(handles, labels, fontsize=legend, loc='lower right', frameon=False, labelcolor='linecolor', markerscale=0, labelspacing = 0.1)

# N against Fe
ax6 = fig.add_subplot(gs[1, 2])
h6 = ax6.hist2d(afe[mask_al], an[mask_al], norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb6 = fig.colorbar(h6[3], ax=ax6, fraction=0.046, pad=0.04)
cb6.set_label('Counts', fontsize=labs)
cb6.ax.tick_params(labelsize=tcks)
bin_centers, bin_x_error_means, bin_y_error_means = bin_mean_errors(feh_g2, nfe_g2, errfeh_g2, errnfe_g2, x_bins)
ax6.errorbar(bin_centers, [1.75] * len(bin_centers), xerr=bin_x_error_means, yerr=bin_y_error_means, 
             ls = 'none', color='blue', label='Mean Errors')
ax6.scatter(feh_al, nfe_al, c='green', alpha=1.0, s=15, label='Dwarf Galaxies')
ax6.scatter(feh_g2[Carina], nfe_g2[Carina], c='red', alpha=0.9, s=50, label='Carina G2s')
ax6.scatter(feh_g2[Draco], nfe_g2[Draco], c='m', alpha=0.9, s=50, label='Draco G2s')
ax6.scatter(feh_g2[Fornax], nfe_g2[Fornax], c='orange', alpha=0.9, s=50, label='Fornax G2s')
ax6.scatter(feh_g2[Sculptor], nfe_g2[Sculptor], c='cyan', alpha=0.9, s=50, label='Sculptor G2s')
ax6.scatter(feh_g2[Ursa_Minor], nfe_g2[Ursa_Minor], c='dodgerblue', alpha=0.9, s=50, label='Ursa Minor G2s')
ax6.set_xlabel('[Fe/H]', size=labs)
ax6.set_xticks(np.arange(-2.5, 1, step=0.5))
ax6.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks, rotation = 45)
ax6.set_ylabel('[N/Fe]', size=labs)
ax6.set_yticks(np.arange(-0.75, 2, step=0.5))
ax6.set_yticklabels(np.arange(-0.75, 2, step=0.5), fontsize=tcks)
ax6.set_box_aspect(1)
#code to remove the marker from the errorbar
handles, labels = ax6.get_legend_handles_labels()
handles = [h[0] if isinstance(h, ErrorbarContainer) else h for h in handles]
ax6.legend(handles, labels, fontsize=legend, loc='upper right', frameon=False, labelcolor='linecolor', markerscale=0, labelspacing = 0.1)

# Mn against Fe
ax7 = fig.add_subplot(gs[2, 0])
h7 = ax7.hist2d(afe[mask_al], amn[mask_al], norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb7 = fig.colorbar(h7[3], ax=ax7, fraction=0.046, pad=0.04)
cb7.set_label('Counts', fontsize=labs)
cb7.ax.tick_params(labelsize=tcks)
bin_centers, bin_x_error_means, bin_y_error_means = bin_mean_errors(feh_g2, mnfe_g2, errfeh_g2, errmnfe_g2, x_bins)
ax7.errorbar(bin_centers, [1] * len(bin_centers), xerr=bin_x_error_means, yerr=bin_y_error_means, 
             ls = 'none', color='blue', label='Mean Errors')
ax7.scatter(feh_al, mnfe_al, c='green', alpha=1.0, s=15, label='Dwarf Galaxies')
ax7.scatter(feh_g2[Carina], mnfe_g2[Carina], c='red', alpha=0.9, s=50, label='Carina G2s')
ax7.scatter(feh_g2[Draco], mnfe_g2[Draco], c='m', alpha=0.9, s=50, label='Draco G2s')
ax7.scatter(feh_g2[Fornax], mnfe_g2[Fornax], c='orange', alpha=0.9, s=50, label='Fornax G2s')
ax7.scatter(feh_g2[Sculptor], mnfe_g2[Sculptor], c='cyan', alpha=0.9, s=50, label='Sculptor G2s')
ax7.scatter(feh_g2[Ursa_Minor], mnfe_g2[Ursa_Minor], c='dodgerblue', alpha=0.9, s=50, label='Ursa Minor G2s')
ax7.set_xlabel('[Fe/H]', size=labs)
ax7.set_xticks(np.arange(-2.5, 1, step=0.5))
ax7.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks, rotation = 45)
ax7.set_ylabel('[Mn/Fe]', size=labs)
ax7.set_yticks(np.arange(-1, 1.5, step=0.5))
ax7.set_yticklabels(np.arange(-1, 1.5, step=0.5), fontsize=tcks)
ax7.set_box_aspect(1)
#code to remove the marker from the errorbar
handles, labels = ax7.get_legend_handles_labels()
handles = [h[0] if isinstance(h, ErrorbarContainer) else h for h in handles]
ax7.legend(handles, labels, fontsize=legend, loc='upper right', frameon=False, labelcolor='linecolor', markerscale=0, labelspacing = 0.1)

# Ca against Fe
#replace nan values from aca with 0
aca_clean = np.nan_to_num(aca[mask_al])

ax8 = fig.add_subplot(gs[2, 1])
h8 = ax8.hist2d(afe[mask_al], aca_clean, norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb8 = fig.colorbar(h8[3], ax=ax8, fraction=0.046, pad=0.04)
cb8.set_label('Counts', fontsize=labs)
cb8.ax.tick_params(labelsize=tcks)
bin_centers, bin_x_error_means, bin_y_error_means = bin_mean_errors(feh_g2, cafe_g2, errfeh_g2, errcafe_g2, x_bins)
ax8.errorbar(bin_centers, [0.8] * len(bin_centers), xerr=bin_x_error_means, yerr=bin_y_error_means, 
             ls = 'none', color='blue', label='Mean Errors')
ax8.scatter(feh_al, cafe_al, c='green', alpha=1.0, s=15, label='Dwarf Galaxies')
ax8.scatter(feh_g2[Carina], cafe_g2[Carina], c='red', alpha=0.9, s=50, label='Carina G2s')
ax8.scatter(feh_g2[Draco], cafe_g2[Draco], c='m', alpha=0.9, s=50, label='Draco G2s')
ax8.scatter(feh_g2[Fornax], cafe_g2[Fornax], c='orange', alpha=0.9, s=50, label='Fornax G2s')
ax8.scatter(feh_g2[Sculptor], cafe_g2[Sculptor], c='cyan', alpha=0.9, s=50, label='Sculptor G2s')
ax8.scatter(feh_g2[Ursa_Minor], cafe_g2[Ursa_Minor], c='dodgerblue', alpha=0.9, s=50, label='Ursa Minor G2s')
ax8.set_xlabel('[Fe/H]', size=labs)
ax8.set_xticks(np.arange(-2.5, 1, step=0.5))
ax8.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks, rotation = 45)
ax8.set_ylabel('[Ca/Fe]', size=labs)
ax8.set_yticks(np.arange(-0.75, 1, step=0.25))
ax8.set_yticklabels(np.arange(-0.75, 1, step=0.25), fontsize=tcks)
ax8.set_box_aspect(1)
#code to remove the marker from the errorbar
handles, labels = ax8.get_legend_handles_labels()
handles = [h[0] if isinstance(h, ErrorbarContainer) else h for h in handles]
ax8.legend(handles, labels, fontsize=legend, loc='upper right', frameon=False, labelcolor='linecolor', markerscale=0, labelspacing = 0.1)

# O against Fe
#replace nan values from aca with 0
ao_clean = np.nan_to_num(ao[mask_al])

ax9 = fig.add_subplot(gs[2, 2])
h9 = ax9.hist2d(afe[mask_al], ao_clean, norm=mpl.colors.LogNorm(), bins=(200, 200), cmap='bone')
cb9 = fig.colorbar(h9[3], ax=ax9, fraction=0.046, pad=0.04)
cb9.set_label('Counts', fontsize=labs)
cb9.ax.tick_params(labelsize=tcks)
bin_centers, bin_x_error_means, bin_y_error_means = bin_mean_errors(feh_g2, ofe_g2, errfeh_g2, errofe_g2, x_bins)
ax9.errorbar(bin_centers, [-0.65] * len(bin_centers), xerr=bin_x_error_means, yerr=bin_y_error_means, 
             ls = 'none', color='blue', label='Mean Errors')
ax9.scatter(feh_al, ofe_al, c='green', alpha=1.0, s=15, label='Dwarf Galaxies')
ax9.scatter(feh_g2[Carina], ofe_g2[Carina], c='red', alpha=0.9, s=50, label='Carina G2s')
ax9.scatter(feh_g2[Draco], ofe_g2[Draco], c='m', alpha=0.9, s=50, label='Draco G2s')
ax9.scatter(feh_g2[Fornax], ofe_g2[Fornax], c='orange', alpha=0.9, s=50, label='Fornax G2s')
ax9.scatter(feh_g2[Sculptor], ofe_g2[Sculptor], c='cyan', alpha=0.9, s=50, label='Sculptor G2s')
ax9.scatter(feh_g2[Ursa_Minor], ofe_g2[Ursa_Minor], c='dodgerblue', alpha=0.9, s=50, label='Ursa Minor G2s')
ax9.set_xlabel('[Fe/H]', size=labs)
ax9.set_xticks(np.arange(-2.5, 1, step=0.5))
ax9.set_xticklabels(np.arange(-2.5, 1, step=0.5), fontsize=tcks, rotation = 45)
ax9.set_ylabel('[O/Fe]', size=labs)
ax9.set_yticks(np.arange(-1, 1.5, step=0.5))
ax9.set_yticklabels(np.arange(-1, 1.5, step=0.5), fontsize=tcks)
ax9.set_box_aspect(1)
#code to remove the marker from the errorbar
handles, labels = ax9.get_legend_handles_labels()
handles = [h[0] if isinstance(h, ErrorbarContainer) else h for h in handles]
ax9.legend(handles, labels, fontsize=legend, loc='lower right', frameon=False, labelcolor='linecolor', markerscale=0, labelspacing = 0.1)

plt.savefig('G2Graphs', bbox_inches='tight')
plt.show()
