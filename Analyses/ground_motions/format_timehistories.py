#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 18:41:54 2023

@author: glavrent
"""

# load packages
import os
import pathlib
import glob
import re  # regular expression package
import pickle
import warnings
# arithmetic libraries
import numpy as np
import numpy.matlib
# statistics libraries
import pandas as pd
# geodesy libraries
import geopy.distance as geopydist
# ploting libraries
import matplotlib as mpl
from matplotlib import pyplot as plt


# %% Define Input Variables
# ======================================
#input/output flag
flag_io = 2
#flag truncate
flag_trnc = True

#gravity
gravity = 9.80665

#component names
cmp_names = ['Z','N','E']

# Input Directories
# ---   ---   ---   ---   ---
#time-histories
if flag_io == 1:  # magnitude 6.5 event
    dir_inp = '../../Raw_files/ground_motions/instrument_corrected/M6.5_0917/ASCII_prc/'
elif flag_io == 2:  # magnitude 6.9 event
    dir_inp = '../../Raw_files/ground_motions/instrument_corrected/M6.9_0918/ASCII_prc/'

#event information
if flag_io == 1:  # magnitude 6.5 event
    eq_info = {'no':111086,'name':'2022 Guanshan',
               'mag':6.5, 'lat':23.08, 'long':121.16,
               'depth':7.3}
elif flag_io == 2:  # magnitude 6.9 event
    eq_info = {'no':111111,'name':'2022 Chihshang',
               'mag':6.9, 'lat':23.14, 'long':121.20,
               'depth':7.0}

#truncation information
trunc_info = {'CWBSN': (25,125), 'EEWS':(20,120), 'SANTA': (5,120), 'TSMIP': (25,135)}

#station information
fn_sta_info     = '../../Raw_files/ground_motions/instrument_corrected/station_info.csv'
fn_sta_metadata = '../../Raw_files/ground_motions/station_metadata.csv'
#station rupture distance
if   flag_io == 1: 
    col_rrup    = 'r_rup'
    col_rjb     = 'r_jb'
    col_rx      = 'r_x'
    col_ry      = 'r_y'
    fn_sta_rrup = '../../Data/ground_motions/dist_metrics/CRF_sta_dist_metrics.csv'
elif flag_io == 2: 
    col_rrup    = 'r_rup'
    col_rjb     = 'r_jb'
    col_rx      = 'r_x'
    col_ry      = 'r_y'
    fn_sta_rrup = '../../Data/ground_motions/dist_metrics/CRF_sta_dist_metrics.csv'

# Output Directories
# ---   ---   ---   ---   ---
if flag_io == 1:
    dir_out = '../../Data/ground_motions/seed_motions/M6.5_0917/'
elif flag_io == 2:
    dir_out = '../../Data/ground_motions/seed_motions/M6.9_0918/'
dir_fig = dir_out + 'figures/'

# %% Processing
# ======================================
#create output directories
if not os.path.isdir(dir_out): pathlib.Path(dir_out).mkdir(parents=True, exist_ok=True)
if not os.path.isdir(dir_fig): pathlib.Path(dir_fig).mkdir(parents=True, exist_ok=True)

#read station information
sta_info = pd.read_csv(fn_sta_info)
sta_meta = pd.read_csv(fn_sta_metadata)
#read rupture distance
sta_rrup = pd.read_csv(fn_sta_rrup)

# find ground-motion files
fnames_gm = np.array(os.listdir(dir_inp))
i_gm = [bool(re.match('.*\.cwb', fn_gm)) for fn_gm in fnames_gm]
fnames_gm = fnames_gm[i_gm]
#number of time histories
n_gm = len(fnames_gm)

#intialize ground-motion info dataframe
df_gm_info = pd.DataFrame(columns=['eq_id','event','mag','hyp_lat','hyp_long','hyp_depth',
                                   'network','station','sta_lat','sta_long', 'sta_elev',
                                   'site_class', 'vs30', 'z1.0', 'kappa', 
                                   'hyp_dist','r_rup','r_jb','r_x','r_y'])

#add pga and file name columns
df_gm_info[['pga_raw_%s'%cmp for cmp in cmp_names]] = np.nan
df_gm_info[['fname_%s'%cmp   for cmp in cmp_names]] = 'nan'

# iterate over ground motions
print('Start Processing Time Histories')
for k, fn_gm in enumerate(fnames_gm):
    print(f'\tProcessing: %14s (%2.i of %i) ...'%(fn_gm,k+1,n_gm))
    #station name
    n_sta = re.search('(.*)\.cwb',fn_gm).group(1)
    #read time history
    #header
    with open(dir_inp + fn_gm) as f:
        contents = f.readlines()
    
    #accelerations
    gm_acc = pd.read_csv(dir_inp + fn_gm, sep="\s+", skiprows=11,
                         header=None, names=['time']+cmp_names)

    #populate ground_motion information
    #event information
    df_gm_info.loc[k,'eq_id'] = eq_info['no']
    df_gm_info.loc[k,'event'] = eq_info['name']
    df_gm_info.loc[k,'mag']   = eq_info['mag']
    df_gm_info.loc[k,'hyp_lat']   = eq_info['lat']
    df_gm_info.loc[k,'hyp_long']  = eq_info['long']
    df_gm_info.loc[k,'hyp_depth'] = eq_info['depth']
    #station information
    s_info = sta_info.loc[sta_info.station==n_sta,:]
    if len(s_info) == 0:
        raise Warning("No station information found.")
    elif len(s_info) > 1:
        raise Warning("Multiple stations with same name.")
    else:
        df_gm_info.loc[k,'network']  = s_info.network.values[0]
        df_gm_info.loc[k,'station']  = s_info.station.values[0]
        df_gm_info.loc[k,'sta_lat']  = s_info.latitude.values[0]
        df_gm_info.loc[k,'sta_long'] = s_info.longitude.values[0]
        df_gm_info.loc[k,'sta_elev'] = s_info.elevation.values[0]/1000 #convert to km
        #compute distance
        d_hyp = geopydist.distance(df_gm_info.loc[k,['hyp_lat','hyp_long']],
                                 df_gm_info.loc[k,['sta_lat','sta_long']]).km
        d_hyp = np.sqrt(d_hyp**2+df_gm_info.loc[k,'hyp_depth']**2)
        df_gm_info.loc[k,'hyp_dist'] = d_hyp 
    #rupture distance information
    s_rup = sta_rrup.loc[sta_rrup.station==n_sta,:]
    if len(s_rup) == 0:
        raise Warning("No station rupture distacne information found.")
    elif len(s_rup) > 1:
        raise Warning("Multiple stations with same name for rupture distacne.")
    else:
        df_gm_info.loc[k,'r_rup'] = s_rup[col_rrup].values[0]
        df_gm_info.loc[k,'r_jb']  = s_rup[col_rjb].values[0]
        df_gm_info.loc[k,'r_x']   = s_rup[col_rx].values[0]
        df_gm_info.loc[k,'r_y']   = s_rup[col_ry].values[0]

    
    #station metadata
    s_meta = sta_meta.loc[sta_meta['Station_code']==n_sta,:]
    if len(s_meta) == 0:
        print("\t\tNo station metadata available.")
    elif len(s_meta) > 1:
        raise Warning("Multiple stations with same name for metadata.")
    else:
        df_gm_info.loc[k,'site_class']  = s_meta['Vs30_classification'].values[0]
        df_gm_info.loc[k,'vs30']        = s_meta['Vs30_sugg'].values[0]        
        df_gm_info.loc[k,'z1.0']        = s_meta['Z1.0_sugg'].values[0]        
        df_gm_info.loc[k,'kappa']       = s_meta['Kappa_sugg'].values[0]        

    #network name
    n_net = df_gm_info.loc[k,'network']

    #extract GM time
    gm_time_raw = contents[2]
    gm_time = [int(e) for e in re.search('^#StartTime.*: (\d+)\/(\d+)\/(\d+)-(\d+):(\d+):(\d+)',gm_time_raw).group(1,2,3,4,5,6)]
    #correct for time zone
    gm_zone = int(re.match('^#StartTime\(GMT(.*)\)', gm_time_raw).group(1)) if bool(re.match('.*\(GMT.*\)', gm_time_raw)) else 0
    gm_time[3] -= gm_zone
    
    #truncate time history    
    if flag_trnc:
        #truncation limits
        trnc_s, trnc_e = trunc_info[n_net]
        #truncation indices
        i_trnc = np.logical_and(gm_acc.time>=trnc_s, gm_acc.time<=trnc_e)
        gm_acc = gm_acc.loc[i_trnc,:]
        #fix time
        gm_acc.time -= trnc_s
        #correct start time
        gm_time[5] += trnc_s
        offset_min = gm_time[5] // 60 #run off minutes
        gm_time[4] += offset_min 
        gm_time[5] -= 60*offset_min 


    #convert gal to g
    for cmp in cmp_names:
        gm_acc.loc[:,cmp] = 0.01 * gm_acc.loc[:,cmp]
        #compute pga
        df_gm_info.loc[k,'pga_raw_'+cmp] = gm_acc.loc[:,cmp].abs().max() / gravity
    
    #save individual components
    for cmp in cmp_names:
        #extract time history component
        gm_acc_c = gm_acc[['time', cmp]]
        #save component time histories
        fn_gm_cmp = '%.4i%.2i%.2i%.2i%.2i%.2i_%s_%s_%s.acc'%tuple(gm_time+[n_net,n_sta,cmp])
        # gm_acc_c.to_csv(dir_out+fn_gm_cmp, index=False, header=False)
        gm_acc_c.to_csv(dir_out+fn_gm_cmp, index=False, header=False, sep=' ', float_format='%8.7f')
        #store filename
        df_gm_info.loc[k,'fname_%s'%cmp] = fn_gm_cmp
        
    #save figures with individual components
    fig, axs = plt.subplots(3, figsize = (40,30)) #create figure
    y_lim = 0.1 * np.ceil(10*np.max([gm_acc.loc[:,cmp].abs().max() for cmp in cmp_names]))
    #iterate over components
    for j, cmp in enumerate(cmp_names):
        #extract time history component
        gm_acc_c = gm_acc[['time', cmp]]
        #plot time histories
        axs[j].plot(gm_acc_c.time, gm_acc_c[cmp], label=cmp, color='black')    
        #fig properties
        axs[j].grid()
        axs[j].legend(loc='upper right', fontsize=60)
        axs[j].set_xlim( (0,axs[j].get_xlim()[1]) )
        #axs[j].set_ylim( np.abs(axs[j].get_ylim()).max() * np.array([-1,1]))
        axs[j].set_ylim( y_lim * np.array([-1,1]) )
        #labels
        axs[j].tick_params(labelsize=50) 
        axs[j].set_ylabel('acc ($m/sec^2$)',    fontsize=60)
        if j == len(cmp_names)-1: axs[j].set_xlabel('Time (sec)', fontsize=60)
    #set title
    axs[0].set_title('Network: '+n_net+', Station: '+n_sta, fontsize=80)
    fig.tight_layout()
    #save time histories
    fn_fig = '%s_%s.png'%(n_net,n_sta)
    fig.savefig( dir_fig+fn_fig, bbox_inches='tight')
    plt.close(fig)

print('Done!')

# %% Output
# ======================================
#earthquake filename
fname_eq = eq_info['name'].replace(' ','_')

#save ground motion info
df_gm_info.to_csv(dir_out+fname_eq+'_gm_info'+'.csv', index=False)

#save ground motion info
df_gm_info4gm_prc = df_gm_info[['event','mag']]
df_gm_info4gm_prc = df_gm_info4gm_prc.join( df_gm_info[['network','station']].agg('_'.join, axis=1).rename('sta') ) 
df_gm_info4gm_prc = df_gm_info4gm_prc.join( df_gm_info[['hyp_dist']] )
#Vs30 info
df_gm_info4gm_prc.loc[range(n_gm),['Vs30']] = -999
#filtering info
df_gm_info4gm_prc.loc[range(n_gm),['h1_HP','h1_LP']] = -999
df_gm_info4gm_prc.loc[range(n_gm),['h2_HP','h2_LP']] = -999
df_gm_info4gm_prc.loc[range(n_gm),['v_HP','v_LP']]   = -999
#time histories info
df_gm_info4gm_prc = df_gm_info4gm_prc.join( df_gm_info[['fname_N','fname_E','fname_Z']] )

#write out file
with open(dir_out+fname_eq+'_gm_info_to_process'+'.csv', 'w') as f:
    f.write('Do not change a column order. Please take a look at Step.0 in R code.\n')
    f.write('Column 1: Earthquake Name,2: Magnitude,3: Station ID,4: Hypocenter Distance (km),5: Vs30 (m/s),6: H1 HP fc,7: H1 LP fc,8: H2 HP fc,9: H2 LP fc,10: V HP fc,11: V LP fc,File Name (Horizontal 1),File Name (Horizontal 2),File Name (Vertical)\n')
df_gm_info4gm_prc.to_csv(dir_out + fname_eq + '_gm_info_to_process' + '.csv',  mode='a', index=False, header=False)
