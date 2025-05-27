#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 00:14:40 2025

@author: glavrent
"""

# load packages
import os
import pathlib
import re
# arithmetic libraries
import numpy as np
#statistics libraries
import pandas as pd
#projection library
import utm
# ploting libraries
import matplotlib as mpl
from matplotlib import pyplot as plt

# %% Define and Load Input Data
# ======================================

#input/output flag
flag_io = 1

#utm zone
zone_number = 51
zone_letter = 'Q'

#filtering info
filter_disp_min = 10.
filter_dist_max = 3000.

#input filenames
fn_smsn_summary           = ['2022_Chihshang_gm_info_summary.csv', '2022_Guanshan_gm_info_summary.csv']
fn_gnss_permanent_summary = ['M6.9_Chihshang_EQ_GNSS_CGS_displacement.csv', 'M6.5_Guanshan_EQ_GNSS_CGS_displacement.csv']
fn_gnss_mobile_summary    = ['M6.9_Chihshang_EQ_GNSS_CGS_mobile_displacement.csv', None]
fn_smsn_gnss              = ['2022_Chihshang_nearby_smsn-gnss.csv', '2022_Guanshan_nearby_smsn-gnss.csv']

#data directory
dir_gis  = '../../../Data/gis/ground_motions/'
dir_gnss = '../../../Raw_files/gnss/'

#output directory
dir_out = '../../../Data/gis/'
dir_fig = dir_out + 'figures/'
    
#read input files
df_smsn = list()
for j in range(len(fn_smsn_summary)):
    #read strong motion datata
    df_smsn.append( pd.read_csv(dir_gis + fn_smsn_summary[j]) )

    #read and merge nearby gnss info
    df_smsn[j] = pd.merge(df_smsn[j], pd.read_csv(dir_gis + fn_smsn_gnss[j]), left_on='station', right_on='SMSN', how='right')
    
    #load gnss data
    df_gnss = pd.read_csv( dir_gnss + fn_gnss_permanent_summary[j] )
    if not fn_gnss_mobile_summary[j] is None:
        df_gnss = pd.concat([df_gnss, pd.read_csv(dir_gnss + fn_gnss_mobile_summary[j])], ignore_index=True)
    df_gnss.rename(columns={'station':'GNSS'}, inplace=True)

    #merge smsn and gnss data
    df_smsn[j] = pd.merge(df_smsn[j], df_gnss, on='GNSS', how='left' )
        
#merge datasets from different events
df_smsn = pd.concat(df_smsn, ignore_index=True)

# %% Processing
# ======================================
#compute utm coordinates for smsn and gnss stations
df_smsn.loc[:,'smsn_utm_x'], df_smsn.loc[:,'smsn_utm_y'], _, _ = utm.from_latlon(df_smsn.sta_lat, 
                                                                                 df_smsn.sta_long, 
                                                                                 force_zone_number=zone_number, 
                                                                                 force_zone_letter=zone_letter)
df_smsn.loc[:,'gnss_utm_x'], df_smsn.loc[:,'gnss_utm_y'], _, _ = utm.from_latlon(df_smsn.latitude, 
                                                                                 df_smsn.longitude, 
                                                                                 force_zone_number=zone_number, 
                                                                                 force_zone_letter=zone_letter)

#compute distance
df_smsn.loc[:,'dist'] = np.linalg.norm(df_smsn.loc[:,['smsn_utm_x','smsn_utm_y']].values - df_smsn.loc[:,['gnss_utm_x','gnss_utm_y']].values, axis=1)


#
df_smsn['tect_slip_azimuth_smsn']  = np.degrees( np.arctan2(df_smsn['tect_slip_E'], df_smsn['tect_slip_N']) ) % 360
df_smsn['tect_slip_azimuth_gnss'] = np.degrees( np.arctan2(df_smsn['long_disp_mean'], df_smsn['lat_disp_mean']) ) % 360


#compute tectonic offset difference
df_smsn['tect_slip_horiz_diff'] = df_smsn['tect_slip_horiz'] - df_smsn['horiz_disp_mean']
df_smsn['tect_slip_vert_diff']  = df_smsn['tect_slip_Z']     - df_smsn['vert_disp_mean']
#compute tectonic offset ratio
df_smsn['tect_slip_horiz_ratio'] = df_smsn['tect_slip_horiz'] / df_smsn['horiz_disp_mean']
df_smsn['tect_slip_vert_ratio']  = df_smsn['tect_slip_Z']     / df_smsn['vert_disp_mean']
#compute tectonic angle difference
df_smsn['tect_slip_azimuth_diff'] = df_smsn.tect_slip_azimuth_gnss - df_smsn.tect_slip_azimuth_smsn 
 
 
# %% Plotting
# ======================================
#create output directories
if not os.path.isdir(dir_out): pathlib.Path(dir_out).mkdir(parents=True, exist_ok=True)
if not os.path.isdir(dir_fig): pathlib.Path(dir_fig).mkdir(parents=True, exist_ok=True)

filter_disp_min = 5
filter_dist_max = 1000.

#comparision horizontal displacement
i_gm = df_smsn.loc[:,'dist'] <= filter_dist_max
fname_fig = 'GNSS-SMSN_horiz_diff'
fig, ax = plt.subplots(figsize = (10,10))
ax.hist(df_smsn.loc[i_gm,'tect_slip_horiz_diff'], bins=np.arange(-25.-1.25, 25.1-1.25, 2.5))
ax.vlines(0, 0, 10, color='k', linewidth=3, zorder=0)
ax.set_xlabel('Horizontal Displacement Difference (cm)',  fontsize=30)
ax.set_ylabel(r'Counts',                                  fontsize=30)
ax.grid(which='both')
ax.set_ylim([0, 10])
ax.set_yticks(range(11))
ax.tick_params(axis='x', labelsize=28)
ax.tick_params(axis='y', labelsize=28)
fig.tight_layout()
fig.savefig( dir_fig + fname_fig + '.png', dpi=300, bbox_inches='tight') 

#comparision vertical displacement
i_gm = df_smsn.loc[:,'dist'] <= filter_dist_max
fname_fig = 'GNSS-SMSN_vert_diff'
fig, ax = plt.subplots(figsize = (9.5,10))
ax.hist(df_smsn.loc[i_gm,'tect_slip_vert_diff'], bins=np.arange(-25.-1.25, 25.1-1.25, 2.5))
ax.vlines(0, 0, 10, color='k', linewidth=3, zorder=0)
ax.set_xlabel('Vertical Displacement Difference (cm)',  fontsize=30)
# ax.set_ylabel(r'Counts',                                  fontsize=30)
ax.grid(which='both')
ax.set_ylim([0, 10])
ax.set_yticks(range(11))
ax.tick_params(axis='x', labelsize=28)
ax.tick_params(axis='y', labelsize=28)
fig.tight_layout()
fig.savefig( dir_fig + fname_fig + '.png', dpi=300, bbox_inches='tight') 

#comparision azimuthal angle
i_gm = np.logical_and(df_smsn.loc[:,'horiz_disp_mean'] >= filter_disp_min, df_smsn.loc[:,'dist'] <= filter_dist_max)
fname_fig = 'GNSS-SMSN_azimuth_diff'
fig, ax = plt.subplots(figsize = (9.5,10))
ax.hist(df_smsn.loc[i_gm,'tect_slip_azimuth_diff'], bins=np.arange(-90., 90.1, 10.))
ax.vlines(0, 0, 10, color='k', linewidth=3, zorder=0)
ax.set_xlabel('Azimuth Angle Difference (deg)', fontsize=30)
# ax.set_ylabel(r'Counts',                        fontsize=30)
ax.grid(which='both')
ax.set_xticks(np.arange(-90, 90.1, 30))
ax.set_xlim([-90, 90])
ax.set_ylim([0, 10])
ax.set_yticks(range(11))
ax.tick_params(axis='x', labelsize=28)
ax.tick_params(axis='y', labelsize=28)
fig.tight_layout()
fig.savefig( dir_fig + fname_fig + '.png', dpi=300, bbox_inches='tight') 
