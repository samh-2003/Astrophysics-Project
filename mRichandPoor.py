# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 10:18:23 2025

@author: Sam
"""

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

#[Fe/H] < -2.2

#define metallicity interval
fmin = -10
fmax = -2.2

        
#Define dwarf galaxies working sample
DG = fits.open('DwarfGalaxies.fits')

ids_al = DG[1].data['APOGEE_ID_1']
teff_al = DG[1].data['TEFF']
logg_al = DG[1].data['LOGG']
mgfe_al = DG[1].data['Mg_Fe']
nfe_al = DG[1].data['N_Fe']
nife_al = DG[1].data['Ni_Fe']
feh_al = DG[1].data['Fe_H']
cafe_al = DG[1].data['Na_Fe']
cfe_al = DG[1].data['C_Fe']
cife_al = DG[1].data['CI_Fe']
ofe_al = DG[1].data['O_Fe']
sife_al = DG[1].data['Si_Fe']
nife_al = DG[1].data['Ni_Fe']
nafe_al = DG[1].data['Na_Fe']
alfe_al = DG[1].data['Al_Fe']
tife_al = DG[1].data['Ti_Fe']
kfe_al = DG[1].data['K_Fe']
vfe_al = DG[1].data['V_Fe']
tife_al = DG[1].data['Ti_Fe']
mnfe_al = DG[1].data['Mn_Fe']
crfe_al = DG[1].data['Cr_Fe']
cufe_al = DG[1].data['Cu_Fe']
cefe_al = DG[1].data['Ce_Fe']
cofe_al = DG[1].data['Co_Fe']
pfe_al = DG[1].data['P_Fe']
sfe_al = DG[1].data['S_Fe']
mgmn_al = DG[1].data['Mg_Fe'] - DG[1].data['Mn_Fe']
simn_al = sife_al - mnfe_al

lz_al = DG[1].data['Lz']
ecc_al = DG[1].data['ecc']
ener_al = DG[1].data['energy']/1.e5
dist_al = DG[1].data['weighted_dist']/1.e3

# Get current size
fig_size = plt.rcParams["figure.figsize"]

# Set figure width to 16 and height to 12
fig_size[0] = 20
fig_size[1] = 20

#Figure settings (font sizes of ticks and labels)
labs=35
tcks=30

#Definition of G2 stars 
lim=0.1
G2 = (alfe_al > lim) 

# These are stars with very high N abundances, absent from GC sample
# Likely debris from w Cen
G3 = nfe_al > 1.5

gs = gridspec.GridSpec(1,2)
gs.update(wspace=0.2, hspace=1) # set the spacing between axes. 

met = ( (feh_al>fmin) & (feh_al < fmax) & (teff_al< 5000) )

plt.subplot(gs[0])
plt.scatter(cfe_al[met],nfe_al[met],c='gray',s=30)
plt.scatter(cfe_al[met&G2],nfe_al[met&G2],c='r',s=50)
plt.scatter(cfe_al[met&G3],nfe_al[met&G3],c='blue',s=50)
#plt.scatter(ct,nt,c='k',alpha=1,s=15)
#plt.scatter(ct[G2],nt[G2],c='r',alpha=1,s=50)
plt.xticks((np.arange(-4,4,step=0.5)),fontsize=tcks)
plt.yticks((np.arange(-4,4,step=0.5)),fontsize=tcks)
plt.ylabel('[N/Fe]',size=labs,labelpad=25)
plt.xlabel('[C/Fe]',size=labs,labelpad=25)
plt.xlim(-1.5,0.8)
plt.ylim(-0.8,2.3)
plt.text(-1.2,-0.5,'Intermediate [Fe/H]',fontsize=25)
plt.tick_params(direction='in',right=True,top=True,length=10)

#plt.plot([-2,2],[lim,lim])
#x = np.arange(-1.5,1,0.1)
#plt.plot(x,a_n*x+b_n, c='blue')



ax = plt.subplot(gs[1])
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
ax.yaxis.set_ticks_position("both")
plt.scatter(mgfe_al[met],alfe_al[met],c='gray',s=30)
plt.scatter(mgfe_al[met&G2],alfe_al[met&G2],c='r',s=50)
plt.scatter(mgfe_al[met&G3],alfe_al[met&G3],c='blue',s=50)
#plt.scatter(mgt,alt,c='k',alpha=1,s=15)
#plt.scatter(mgt[G2],alt[G2],c='r',alpha=1,s=25)
plt.xticks((np.arange(-4,4,step=0.2)),fontsize=tcks)
plt.yticks((np.arange(-4,4,step=0.5)),fontsize=tcks)
plt.ylabel('[Al/Fe]',size=labs,labelpad=25)
plt.xlabel('[Mg/Fe]',size=labs,labelpad=25)
plt.xlim(-0.5,0.7)
plt.ylim(-1.0,2.0)
plt.tick_params(direction='in',right=True,top=True,length=10,labelright=True,labelleft=False)
#plt.plot([-2,2],[lim,lim])
#x = np.arange(-1.5,1,0.1)
#plt.plot(x,a_al*x+b_al, c='blue')




#plt.savefig('train_'+ pt_name[ind] + '.png',format='png',dpi=300,bbox_inches='tight')