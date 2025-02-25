# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:27:29 2025

@author: samhe
"""

#### -1.3 < [Fe/H] < -1.0


#import relevant modules
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

#Import GC VAC file
g = fits.open('GC_members_VAC-v1_1.fits')
gc_name = g[1].data['GC_NAME']
gc_id = g[1].data['APOGEE_ID']
gc_tid = g[1].data['TARGET_ID']
gc_gid = g[1].data['GAIAEDR3_SOURCE_ID']
gc_ra = g[1].data['RA']
gc_dec = g[1].data['DEC']
gc_glon = g[1].data['GLON']
gc_glat = g[1].data['GLAT']
gc_vscat = g[1].data['VSCATTER']
gc_nvis = g[1].data['NVISITS']
gc_teff = g[1].data['TEFF']
gc_teff_e = g[1].data['TEFF_ERR']
gc_logg = g[1].data['LOGG']
gc_logg_e = g[1].data['LOGG_ERR']
gc_sn = g[1].data['SNR']
gc_rv = g[1].data['VHELIO_AVG']


gc_feh = g[1].data['Fe_H']
gc_mgfe = g[1].data['Mg_Fe']
gc_cfe = g[1].data['C_Fe']
gc_cife = g[1].data['CI_Fe']
gc_ofe = g[1].data['O_Fe']
gc_nfe = g[1].data['N_Fe']
gc_cafe = g[1].data['Ca_Fe']
gc_sife = g[1].data['Si_Fe']
gc_nife = g[1].data['Ni_Fe']
gc_alfe = g[1].data['Al_Fe']
gc_tife = g[1].data['Ti_Fe']
gc_cofe = g[1].data['Co_Fe']
gc_sfe = g[1].data['S_Fe']
gc_kfe = g[1].data['K_Fe']
gc_pfe = g[1].data['P_Fe']
gc_cefe = g[1].data['Ce_Fe']
gc_vfe = g[1].data['V_Fe']
gc_crfe = g[1].data['Cr_Fe']
gc_nafe = g[1].data['Na_Fe']
gc_mnfe = g[1].data['Mn_Fe']
gc_cufe = g[1].data['Cu_Fe']
gc_sflag = g[1].data['STARFLAG']
gc_aflag = g[1].data['ASPCAPFLAG']
gc_saflag = g[1].data['ASPCAPFLAGS']
gc_ealfe = g[1].data['Al_FE_ERR']
gc_emgfe = g[1].data['MG_FE_ERR']
gc_emnfe = g[1].data['MN_FE_ERR']
gc_enfe = g[1].data['N_FE_ERR']
gc_enife = g[1].data['NI_FE_ERR']
gc_esife = g[1].data['SI_FE_ERR']
gc_ecfe = g[1].data['C_FE_ERR']

#subtract duplicates using pandas

glisttup=list(zip(gc_id, gc_tid, gc_gid, gc_ra, gc_dec, gc_glon, gc_glat, gc_vscat, gc_nvis, \
            gc_teff, gc_teff_e, gc_logg, gc_logg_e, gc_sn, gc_rv, gc_feh, gc_mgfe, \
            gc_cfe, gc_cife, gc_ofe, gc_nfe, gc_cafe, gc_sife, gc_nife, gc_alfe, gc_tife, \
            gc_cofe, gc_sfe, gc_kfe, gc_pfe, gc_cefe, gc_vfe, gc_crfe, gc_nafe, gc_mnfe, \
            gc_cufe, gc_sflag, gc_aflag, gc_ecfe, gc_enfe, gc_enife, gc_esife, gc_emnfe, \
            gc_emgfe, gc_ealfe, gc_name, gc_saflag))

g_df=pd.DataFrame(glisttup, \
        columns=['gc_id', 'gc_tid', 'gc_gid', 'gc_ra', 'gc_dec', 'gc_glon', 'gc_glat', \
        'gc_vscat', 'gc_nvis', 'gc_teff', 'gc_teff_e', 'gc_logg', 'gc_logg_e', \
		'gc_sn', 'gc_rv', 'gc_feh', 'gc_mgfe', 'gc_cfe', 'gc_cife', 'gc_ofe', 'gc_nfe', \
        'gc_cafe', 'gc_sife', 'gc_nife', 'gc_alfe', 'gc_tife', 'gc_cofe', 'gc_sfe', 'gc_kfe', \
        'gc_pfe', 'gc_cefe', 'gc_vfe', 'gc_crfe', 'gc_nafe', 'gc_mnfe', 'gc_cufe', 'gc_sflag', \
		'gc_aflag', 'gc_ecfe', 'gc_enfe', 'gc_enife', 'gc_esife', 'gc_emnfe', 'gc_emgfe', \
        'gc_ealfe', 'gc_name', 'gc_saflag'])

#remove duplicates, keeping the records with the highest SNR value
print('Number of stars before removing duplicates: '+str(len(g_df)))
g_df=g_df.sort_values(by=['gc_id','gc_sn'],ascending=False)
g_df=g_df[~g_df['gc_id'].duplicated(keep='first')]
print('Number of stars after removing duplicates: '+str(len(g_df)))

g_df=g_df.reset_index(drop=True)

#convert back to numpy
gc_name = g_df['gc_name'].to_numpy()
gc_id = g_df['gc_id'].to_numpy()
gc_ra = g_df['gc_ra'].to_numpy()
gc_dec = g_df['gc_dec'].to_numpy()
gc_glon = g_df['gc_glon'].to_numpy()
gc_glat = g_df['gc_glat'].to_numpy()
gc_vscat = g_df['gc_vscat'].to_numpy()
gc_nvis = g_df['gc_nvis'].to_numpy()
gc_teff = g_df['gc_teff'].to_numpy()
gc_logg = g_df['gc_logg'].to_numpy()
gc_sn = g_df['gc_sn'].to_numpy()
gc_rv = g_df['gc_rv'].to_numpy()
gc_feh = g_df['gc_feh'].to_numpy()
gc_mgfe = g_df['gc_mgfe'].to_numpy()
gc_cfe = g_df['gc_cfe'].to_numpy()
gc_cife = g_df['gc_cife'].to_numpy()
gc_ofe = g_df['gc_ofe'].to_numpy()
gc_nfe = g_df['gc_nfe'].to_numpy()
gc_cafe = g_df['gc_cafe'].to_numpy()
gc_sife = g_df['gc_sife'].to_numpy()
gc_nife = g_df['gc_nife'].to_numpy()
gc_alfe = g_df['gc_alfe'].to_numpy()
gc_tife = g_df['gc_tife'].to_numpy()
gc_cofe = g_df['gc_cofe'].to_numpy()
gc_sfe = g_df['gc_sfe'].to_numpy()
gc_kfe = g_df['gc_kfe'].to_numpy()
gc_pfe = g_df['gc_pfe'].to_numpy()
gc_cefe = g_df['gc_cefe'].to_numpy()
gc_vfe = g_df['gc_vfe'].to_numpy()
gc_crfe = g_df['gc_crfe'].to_numpy()
gc_nafe = g_df['gc_nafe'].to_numpy()
gc_mnfe = g_df['gc_mnfe'].to_numpy()
gc_cufe = g_df['gc_cufe'].to_numpy()
gc_sflag = g_df['gc_sflag'].to_numpy()
gc_aflag = g_df['gc_aflag'].to_numpy()
gc_saflag = g_df['gc_saflag'].to_numpy()

gc_ecfe = g_df['gc_ecfe'].to_numpy()  
gc_enfe = g_df['gc_enfe'].to_numpy()
gc_enife = g_df['gc_enife'].to_numpy()
gc_esife = g_df['gc_esife'].to_numpy()
gc_emnfe = g_df['gc_emnfe'].to_numpy()  
gc_emgfe = g_df['gc_emgfe'].to_numpy()  
gc_ealfe = g_df['gc_ealfe'].to_numpy()


#Identify and remove NaN's
mask_nan = (np.isfinite(gc_alfe) & np.isfinite(gc_feh) & np.isfinite(gc_nfe) \
            & np.isfinite(gc_cfe) & np.isfinite(gc_mgfe))
mask_g = ( (gc_sn >= 50) & (gc_teff > 3000) & (gc_teff < 5000) & \
        (gc_logg < 3.6) & (gc_logg > -1) )
        
#Define dwarf galaxies working sample
DG = fits.open('DwarfGalaxies.fits')

ids_al = DG[1].data['APOGEE_ID_1']
afield_al = DG[1].data['FIELD']
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
a_n = -0.5
b_n = 0.3

lim_n = a_n*cfe_al + b_n

G2 = (nfe_al > lim_n) 

lim_n = a_n*gc_cfe + b_n

gc_G2 = (gc_nfe > lim_n)  


gs = gridspec.GridSpec(2,2)
gs.update(wspace=0.2, hspace=0.2) # set the spacing between axes. 

#define metallicity interval
fmin = -1.3
fmax = -1.0

#Define metallicity ranges
met = ( (feh_al>fmin) & (feh_al < fmax) & (teff_al< 5000) )
gc_met = ( (gc_feh>fmin) & (gc_feh < fmax) & (gc_teff < 5000) )

plt.subplot(gs[0])
plt.rc('text', usetex=True)
plt.suptitle(r'$\underline{-1.0 > [Fe/H] > -1.3}$',y = 0.94, fontsize = 50, usetex = True)
plt.scatter(gc_cfe[gc_met&mask_nan&mask_g],gc_nfe[gc_met&mask_nan&mask_g],c='k',s=30,alpha=0.7, label = 'GC VAC')
#plt.scatter(cfe_al[met&G3],nfe_al[met&G3],c='blue',s=60)
plt.scatter(gc_cfe[gc_met&mask_nan&mask_g&gc_G2],gc_nfe[gc_met&mask_nan&mask_g&gc_G2],c='r',alpha=0.7,s=60, label = 'GC VAC G2', marker = 'x')
#plt.scatter(gc_cfe,gc_nfe,c='k',alpha=1,s=15)
plt.scatter(cfe_al[met&G2],nfe_al[met&G2],c='g',s=60, label = 'Dwarf Galaxy G2')
plt.xticks((np.arange(-1.5,1.5,step=0.5)),fontsize=tcks)
plt.yticks((np.arange(-1,2.5,step=0.5)),fontsize=tcks)
plt.ylabel('[N/Fe]',size=labs,labelpad=25)
plt.xlabel('[C/Fe]',size=labs,labelpad=25)
plt.xlim(-1.5,0.8)
plt.ylim(-0.8,2.3)
plt.gca().set_box_aspect(1)
#plt.text(-1.2,-0.65,'Intermediate [Fe/H]',fontsize=25)
plt.legend(fontsize = 15, loc='upper right')
plt.tick_params(direction='in',right=True,top=True,length=10)
#plt.plot([-1.5, -1.2],[0, 0.6], c='blue')



ax = plt.subplot(gs[1])
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
ax.yaxis.set_ticks_position("both")
plt.scatter(gc_mgfe[gc_met&mask_nan&mask_g],gc_alfe[gc_met&mask_nan&mask_g],c='k',s=30,alpha=0.7, label = 'GC VAC')
#plt.scatter(mgfe_al[met&G3],alfe_al[met&G3],c='blue',s=60)
plt.scatter(gc_mgfe[gc_met&mask_nan&mask_g&gc_G2],gc_alfe[gc_met&mask_nan&mask_g&gc_G2],c='r',alpha=0.7,s=60, label = 'GC VAC G2', marker = 'x')
#plt.scatter(gc_mgfe,gc_alfe,c='k',alpha=1,s=15)
plt.scatter(mgfe_al[met&G2],alfe_al[met&G2],c='g',s=60, label = 'Dwarf Galaxy G2')
plt.xticks((np.arange(-4,4,step=0.2)),fontsize=tcks)
plt.yticks((np.arange(-4,4,step=0.5)),fontsize=tcks)
plt.ylabel('[Al/Fe]',size=labs,labelpad=25)
plt.xlabel('[Mg/Fe]',size=labs,labelpad=25)
plt.xlim(-0.5,0.7)
plt.ylim(-1.0,2.0)
plt.gca().set_box_aspect(1)
plt.legend(fontsize = 15, loc='upper right')
plt.tick_params(direction='in',right=True,top=True,length=10,labelright=True,labelleft=False)




plt.subplot(gs[2])
plt.scatter(gc_alfe[gc_met&mask_nan&mask_g],gc_nfe[gc_met&mask_nan&mask_g],c='k',s=30,alpha=0.7, label = 'GC VAC')
#plt.scatter(cfe_al[met&G3],nfe_al[met&G3],c='blue',s=60)
plt.scatter(gc_alfe[gc_met&mask_nan&mask_g&gc_G2],gc_nfe[gc_met&mask_nan&mask_g&gc_G2],c='r',alpha=0.7,s=60, label = 'GC VAC G2', marker = 'x')
#plt.scatter(gc_cfe,gc_nfe,c='k',alpha=1,s=15)
plt.scatter(alfe_al[met&G2],nfe_al[met&G2],c='g',s=60, label = 'Dwarf Galaxy G2')
plt.xticks((np.arange(-1.5,1.5,step=0.5)),fontsize=tcks)
plt.yticks((np.arange(-1,2.5,step=0.5)),fontsize=tcks)
plt.ylabel('[N/Fe]',size=labs,labelpad=25)
plt.xlabel('[Al/Fe]',size=labs,labelpad=25)
plt.xlim(-1.5,0.8)
plt.ylim(-0.8,2.3)
plt.gca().set_box_aspect(1)
#plt.text(-1.2,-0.65,'Intermediate [Fe/H]',fontsize=25)
plt.legend(fontsize = 15, loc='upper right')
plt.tick_params(direction='in',right=True,top=True,length=10)




ax = plt.subplot(gs[3])
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
ax.yaxis.set_ticks_position("both")
plt.scatter(gc_mgfe[gc_met&mask_nan&mask_g],gc_nfe[gc_met&mask_nan&mask_g],c='k',s=30,alpha=0.7, label = 'GC VAC')
#plt.scatter(mgfe_al[met&G3],alfe_al[met&G3],c='blue',s=60)
plt.scatter(gc_mgfe[gc_met&mask_nan&mask_g&gc_G2],gc_nfe[gc_met&mask_nan&mask_g&gc_G2],c='r',alpha=0.7,s=60, label = 'GC VAC G2', marker = 'x')
#plt.scatter(gc_mgfe,gc_alfe,c='k',alpha=1,s=15)
plt.scatter(mgfe_al[met&G2],cfe_al[met&G2],c='g',s=60, label = 'Dwarf Galaxy G2')
plt.xticks((np.arange(-4,4,step=0.2)),fontsize=tcks)
plt.yticks((np.arange(-4,4,step=0.5)),fontsize=tcks)
plt.ylabel('[N/Fe]',size=labs,labelpad=25)
plt.xlabel('[Mg/Fe]',size=labs,labelpad=25)
plt.xlim(-0.5,0.7)
plt.ylim(-1.0,2.0)
plt.gca().set_box_aspect(1)
plt.legend(fontsize = 15, loc='upper right')
plt.tick_params(direction='in',right=True,top=True,length=10,labelright=True,labelleft=False)


print(afield_al[met&G2])
#print(afield_al[met])


#plt.savefig('train_'+ pt_name[ind] + '.png',format='png',dpi=300,bbox_inches='tight')

#Save newly found G2 stars to G2 catalogue and check if its done correctly
#define list sizes
#Table.
Table.max_lines = 50
#Table.
Table.max_width = 100


#define the hdu accessed and selected data
hdu = DG[1]
from astropy.table import vstack
selected_data = Table(hdu.data[met&G2])
hdul = fits.open('DwarfG2s.fits')
target_data = Table(hdul[1].data)
hdul.close()


#check table before appending
print('Table Before Appending:')
hdul = fits.open('DwarfG2s.fits')
for hdu in hdul:
    if isinstance(hdu, fits.BinTableHDU):
        table = Table(hdu.data)
        print(table)
        break
hdul.close()

#stack both data sets and write them to the fits file
combined_data = vstack([target_data, selected_data])
combined_data.write('DwarfG2s.fits', overwrite = True)

#check new values have been appended
print('Table After Appending:')
hdul = fits.open('DwarfG2s.fits')
for hdu in hdul:
    if isinstance(hdu, fits.BinTableHDU):
        table = Table(hdu.data)
        print(table)
        break
hdul.close()