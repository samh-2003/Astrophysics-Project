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

hdu = G2[1]
ids_g2 = G2[1].data['APOGEE_ID_1']
afield_g2 = G2[1].data['FIELD']
mgfe_g2 = G2[1].data['Mg_Fe']
nfe_g2 = G2[1].data['N_Fe']
cfe_g2 = G2[1].data['C_Fe']
alfe_g2 = G2[1].data['Al_Fe']

#Open GC field
gc = fits.open('GC_members_VAC-v1_1.fits')
mgfe_gc = gc[1].data['Mg_Fe']
nfe_gc = gc[1].data['N_Fe']
cfe_gc = gc[1].data['C_Fe']
alfe_gc = gc[1].data['Al_Fe']

#create mask to find the high carbon 'G2s'
#carbonBad = (cfe_g2 > -0.25)

#save carbon bad stars to a seperate file
#selected_data = hdu.data[carbonBad]
#image_hdu = fits.BinTableHDU(data=selected_data, name='Dwarf Galaxy Data')
#print(image_hdu)
#image_hdu.writeto('BadCarbon.fits', overwrite=True)

#import bad carbon stars
BC = fits.open('BadCarbon.fits')
mgfe_bc = BC[1].data['Mg_Fe']
nfe_bc = BC[1].data['N_Fe']
cfe_bc = BC[1].data['C_Fe']
alfe_bc = BC[1].data['Al_Fe']

# Check for NaN values in mgfe_bc and alfe_bc
nan_mgfe_bc = np.isnan(mgfe_bc)
nan_alfe_bc = np.isnan(alfe_bc)

# Count the number of NaN values
num_nan_mgfe_bc = np.sum(nan_mgfe_bc)
num_nan_alfe_bc = np.sum(nan_alfe_bc)

# Print the results
print(f"Number of NaN values in mgfe_bc: {num_nan_mgfe_bc}")
print(f"Number of NaN values in alfe_bc: {num_nan_alfe_bc}")

# Check for NaN values in the G2 sample overall
nan_mgfe_g2 = np.isnan(mgfe_g2)
nan_alfe_g2 = np.isnan(alfe_g2)

# Count the number of NaN values
num_nan_mgfe_g2 = np.sum(nan_mgfe_g2)
num_nan_alfe_g2 = np.sum(nan_alfe_g2)

# Print the results
print(f"Number of NaN values in mgfe_g2: {num_nan_mgfe_g2}")
print(f"Number of NaN values in alfe_g2: {num_nan_alfe_g2}")


#Graphs
fig = plt.figure(figsize=(30, 30))
gs = gridspec.GridSpec(1, 2, figure=fig)

plt.rc('text', usetex=True)
plt.suptitle(r'$\underline{Possible\;Incorrect\;G2\;Identifications}$', fontsize = 50, y = 0.725)

# C against N
ax1 = fig.add_subplot(gs[0, 0])
ax1.scatter(cfe_gc, nfe_gc, c='black', alpha=1.0, s=15, label='GC field')
ax1.scatter(cfe_g2, nfe_g2, c='green', alpha=0.9, s=50, label='Other incorrect G2s')
ax1.scatter(cfe_bc, nfe_bc, c='red', alpha=0.9, s=50, label='Possible incorrect G2s')
ax1.set_xlabel('[C/Fe]', size=labs)
ax1.set_xticks(np.arange(-1.5, 1.5, step=0.5))
ax1.set_xticklabels(np.arange(-1.5, 1.5, step=0.5), fontsize=tcks, rotation = 45)
ax1.set_ylabel('[N/Fe]', size=labs)
ax1.set_yticks(np.arange(-1, 2.5, step=0.5))
ax1.set_yticklabels(np.arange(-1, 2.5, step=0.5), fontsize=tcks)
ax1.set_box_aspect(1)
ax1.legend(fontsize=30, loc='upper right')

# Mg against Al
ax2 = fig.add_subplot(gs[0, 1])
ax2.scatter(mgfe_gc, alfe_gc, c='black', alpha=1.0, s=15, label='GC field')
ax2.scatter(mgfe_g2, alfe_g2, c='green', alpha=0.9, s=50, label='Other G2s')
ax2.scatter(mgfe_bc, alfe_bc, c='red', alpha=0.9, s=50, label='Possible incorrect G2s')
ax2.set_xlabel('[Mg/Fe]', size=labs)
ax2.set_xticks(np.arange(-1, 1, step=0.25))
ax2.set_xticklabels(np.arange(-1, 1, step=0.25), fontsize=tcks, rotation = 45)
ax2.set_ylabel('[Al/Fe]', size=labs)
ax2.set_yticks(np.arange(-2, 2.5, step=0.5))
ax2.set_yticklabels(np.arange(-2, 2.5, step=0.5), fontsize=tcks)
ax2.set_box_aspect(1)
ax2.legend(fontsize=30, loc='lower right')

plt.show()