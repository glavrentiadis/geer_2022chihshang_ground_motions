#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 19:09:44 2023

@author: glavrent
"""

# load packages
import os
import pathlib
# arithmetic libraries
import numpy as np
#statistics libraries
import pandas as pd
# ploting libraries
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#user libraries
exec(open('../../python_lib/plotting/pylib_contour_plots.py').read())

# %% Define and Load Input Data
# ======================================
#dataset flag
# flag_ds = 1: entire dataset
# flag_ds = 2: processed dataset
flag_ds = 2

#taiwan flatfile
if flag_ds == 1:   #input dataset (entire dataset)
    fname_M65_Guanshan  = '../../../Data/ground_motions/seed_motions/M6.5_0917/2022_Guanshan_gm_info.csv'
    fname_M69_Chihshang = '../../../Data/ground_motions/seed_motions/M6.9_0918/2022_Chihshang_gm_info.csv'
elif flag_ds == 2: #input dataset (processed dataset)
    fname_M65_Guanshan  = '../../../Data/ground_motions/seed_motions/2022_Guanshan_gm_info_for_processing.csv'
    fname_M69_Chihshang = '../../../Data/ground_motions/seed_motions/2022_Chihshang_gm_info_for_processing.csv'

#finte fault trace
fname_LVF_trace = '../../../Raw_files/gis/LVF_finite_fault_trace.csv'
fname_LVF_box   = '../../../Raw_files/gis/LVF_finite_fault_boundary.csv'
fname_CRF_trace = '../../../Raw_files/gis/CRF_finite_fault_trace.csv'
fname_CRF_box   = '../../../Raw_files/gis/CRF_finite_fault_boundary.csv'

#NGA2 flatfile
fname_flatfile_NGA2ASK14 = '../../../Raw_files/nga_data/resid_T0.200.out2.txt'
fname_flatfile_NGA2coor  = '../../../Raw_files/nga_data/Updated_NGA_West2_Flatfile_coordinates.csv'
#NGA3 flatfile
fname_flatfile_NGA3 = '../../../Raw_files/nga_data/NGAW3_FAS_Data_F_20241116_karen_v7.csv'

#output directory
dir_out = '../../../Data/ground_motions/summary_general/'
dir_fig = dir_out + 'figures/'

#read dataframes
df_gm_guanshan  = pd.read_csv(fname_M65_Guanshan)
df_gm_chihshang = pd.read_csv(fname_M69_Chihshang)

#read fualt traces
df_flt_LVF_trace = pd.read_csv(fname_LVF_trace)
df_flt_CRF_trace = pd.read_csv(fname_CRF_trace)
df_flt_LVF_box   = pd.read_csv(fname_LVF_box)
df_flt_CRF_box   = pd.read_csv(fname_CRF_box)

#read NGAWest2
df_flatfile_NGA2ASK14 = pd.read_csv(fname_flatfile_NGA2ASK14, delim_whitespace=True)
df_flatfile_NGA2coor  = pd.read_csv(fname_flatfile_NGA2coor)
df_flatfile_NGA2 = pd.merge(df_flatfile_NGA2ASK14, df_flatfile_NGA2coor, left_on='recID', right_on='Record Sequence Number')
del df_flatfile_NGA2ASK14, df_flatfile_NGA2coor

#read NGAWest3
df_flatfile_NGA3 = pd.read_csv(fname_flatfile_NGA3)
df_flatfile_NGA3 = df_flatfile_NGA3.loc[~np.isin(df_flatfile_NGA3.event_id, [3591, 3592]),:].reset_index(drop=True)

# %% Summary
# ======================================
msg = list()
msg.append(f'Summary compiled stong-motions:')
msg.append(f'\tTotal number of records: %i'%(len(df_gm_guanshan)+len(df_gm_chihshang)))
msg.append(f'\tM6.5 Guanshan records: %i'%(len(df_gm_guanshan)))
msg.append(f'\tM6.9 Chihshang records: %i'%(len(df_gm_chihshang)))
msg.append(f'\tR_rup < 10 km records: %i'%( (df_gm_guanshan.r_rup <= 10).sum() + (df_gm_chihshang.r_rup <= 10).sum()))
msg.append(f'\tR_hyp < 10 km records: %i'%( (df_gm_guanshan.hyp_dist <= 10).sum() + (df_gm_chihshang.hyp_dist <= 10).sum()))
msg.append(f'Summary NGAW2:')
msg.append(f'\tR_rup < 10 km records: %i'%( (df_flatfile_NGA2.loc[df_flatfile_NGA2.mag>=6.5,'Rrup'] <= 10).sum()))
msg.append(f'Summary NGAW3:')
msg.append(f'\tR_rup < 10 km records: %i'%( (df_flatfile_NGA3.loc[df_flatfile_NGA3.magnitude>=6.5,'rrup_m'] <= 10).sum()))

#print summary
for m in msg:
    print(m)

#save summary
fname_out = 'summary_gm_taiwan' + '.txt'
with open(dir_out+fname_out, 'w') as f:
    for m in msg:
        f.write(m+'\n')


# %% Plotting
# ======================================
#create output directories
if not os.path.isdir(dir_out): pathlib.Path(dir_out).mkdir(parents=True, exist_ok=True)
if not os.path.isdir(dir_fig): pathlib.Path(dir_fig).mkdir(parents=True, exist_ok=True)

#magnitude-distance distribution
fname_fig = 'distribution_M-R_taiwan'
fig, ax = plt.subplots(figsize = (10,10))
ax.plot(df_gm_chihshang.r_rup,    df_gm_chihshang.mag,  's', markersize=12,               label='2021 M6.9 Chihshang')
ax.plot(df_gm_guanshan.r_rup,     df_gm_guanshan.mag,   'd', markersize=12,               label='2021 M6.5 Guanshan')
# ax.plot(df_gm_chihshang.hyp_dist, df_gm_chihshang.mag,  's', markersize=12,               label='2021 M6.9 Chihshang')
# ax.plot(df_gm_guanshan.hyp_dist,  df_gm_guanshan.mag,   'd', markersize=12,               label='2021 M6.5 Guanshan')
#edit figure
ax.set_xlim([0., 50.])
ax.set_ylim([3., 8.])
ax.set_xlabel(r'$R_{rup}$ ($km$)', fontsize=35)
# ax.set_xlabel(r'Hypocenter Distance (km)',  fontsize=35)
ax.set_ylabel(r'Magnitude',                 fontsize=35)
ax.legend(loc='lower right', fontsize=35)
ax.grid(which='both')
ax.tick_params(axis='x', labelsize=30)
ax.tick_params(axis='y', labelsize=30)
fig.tight_layout()
fig.savefig(dir_fig + fname_fig + '.png')

#magnitude-distance distribution
fname_fig = 'distribution_M-R_taiwan-NGAW2-NGA3'
fig, ax = plt.subplots(figsize = (10,10))
ax.plot(df_flatfile_NGA2.Rrup,    df_flatfile_NGA2.mag,       'o', markersize=5,  color='black', label='NGAWest2', zorder=3)
ax.plot(df_flatfile_NGA3.rrup_m,  df_flatfile_NGA3.magnitude, 'o', markersize=3,  color='gray',  label='NGAWest3', zorder=2)
ax.plot(df_gm_chihshang.r_rup,    df_gm_chihshang.mag,        's', markersize=12,                label='2021 M6.9 Chishang', zorder=4)
ax.plot(df_gm_guanshan.r_rup,     df_gm_guanshan.mag,         'd', markersize=12,                label='2021 M6.5 Guanshan', zorder=4)
# ax.plot(df_gm_chihshang.hyp_dist, df_gm_chihshang.mag,      's', markersize=12,                label='2021 M6.9 Chihshang')
# ax.plot(df_gm_guanshan.hyp_dist,  df_gm_guanshan.mag,       'd', markersize=12,                label='2021 M6.5 Guanshan')
#edit figure
ax.set_xscale('log')
ax.set_xlim([0., 500.])
ax.set_ylim([2., 8.])
ax.set_xlabel(r'$R_{rup}$ ($km$)', fontsize=35)
ax.set_ylabel(r'Magnitude',        fontsize=35)
ax.legend(loc='lower left', fontsize=35)
ax.grid(which='both')
ax.tick_params(axis='x', labelsize=30)
ax.tick_params(axis='y', labelsize=30)
fig.tight_layout()
fig.savefig(dir_fig + fname_fig + '.png')


#eq and sta locations
fname_fig = 'eq_sta_locations'
fig, ax, data_crs, gl = PlotMap(flag_grid=True)
#plot earthquake and station locations
ax.plot(df_gm_chihshang['hyp_long'].values, df_gm_chihshang['hyp_lat'].values, '*', 
        transform=data_crs, markersize=35, zorder=12, label='2021 M6.9 Chishang')
ax.plot(df_gm_guanshan['hyp_long'].values,  df_gm_guanshan['hyp_lat'].values, '*', 
        transform=data_crs, markersize=35, zorder=12, label='2021 M6.5 Guanshan')
#plot fault traces
ax.plot(df_flt_CRF_trace.longitude, df_flt_CRF_trace.latitude, '-k',  linewidth=3., label='Central Valley Fault')
ax.plot(df_flt_CRF_box.longitude, df_flt_CRF_box.latitude,     ':k', linewidth=1.)
ax.plot(df_flt_LVF_trace.longitude, df_flt_LVF_trace.latitude, '--k',  linewidth=3., label='Long Valley Fault')
ax.plot(df_flt_LVF_box.longitude, df_flt_LVF_box.latitude,     ':k', linewidth=1.)
#plot stations
ax.plot(df_gm_guanshan['sta_long'].values,   df_gm_guanshan['sta_lat'].values, 'o', 
        transform=data_crs, markersize=8, zorder=3, color='black', label='Strong Motion Stations')
#edit figure properties
gl.xlabel_style = {'size': 25}
gl.ylabel_style = {'size': 25}
gl.top_labels   = False
gl.right_labels = False 
# ax.legend(fontsize=25, loc='upper left')
# ax.legend(fontsize=25, loc='lower right')
ax.legend(fontsize=25, loc='upper left', bbox_to_anchor=(1, 0.5))
gl.xlocator = mticker.FixedLocator([121.0, 121.2, 121.4, 121.6, 121.8])
gl.ylocator = mticker.FixedLocator([22.6, 23.0, 23.4, 23.8])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
# ax.set_xticks([121.0, 121.4, 121.8])
# ax.set_xticklabels([])
ax.set_xlim([120.8, 122.0])
ax.set_ylim([22.5,  24.0])
# ax.set_xlim([120.8, 121.8])
# ax.set_ylim([22.5,  24.0])
#save figure
fig.tight_layout()
fig.savefig( dir_fig + fname_fig + '.png', dpi=300, bbox_inches='tight') 
