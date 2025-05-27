#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 22:19:24 2023

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
# ploting libraries
import matplotlib as mpl
from matplotlib import pyplot as plt

# %% Define and Load Input Data
# ======================================
#flag components
flag_rot = 0

#gravity
gravity = 9.80665

#components
cmp_names = ['Z','FN','FP'] if flag_rot else ['Z','N','E']
cmp_names = ['N','E','Z']

#ground motion data
dir_gm = '../../../Data/ground_motions/vel_pulse/' + 'M6.9_0918/'
n_sta = ['HWA054','HWA073']

#axis limits
ylim_acc = [7.5, 7.5]
ylim_vel = [100., 100.]
ylim_dis = [100., 120.]
#line parameters
lw_map = [1.5,2]
# lc_map = ['#ff7f0e','#1f77b4']
lc_map = ['black','#1f77b4']
ls_map  = ['-','-']

#label information
label = {'processed':'Time History',
         'vel_pulse':'Velocity Pulse'}

#filename processing info
fn_prc_info = '2022_Chihshang_gm_info_vel_pulse.csv'

#figure directory
dir_out = '../../../Data/ground_motions/vel_pulse/summary/'
dir_fig = dir_out + 'figures/'

# %% Plotting
# ======================================
#create output directories
if not os.path.isdir(dir_fig): pathlib.Path(dir_fig).mkdir(parents=True, exist_ok=True)

#processing info
df_prc_info = pd.read_csv(dir_gm+fn_prc_info)

#iterate over selected ground motions
for k, n_s in enumerate(n_sta):
    #read time histories
    # - - - - - - - - - - - - 
    #ground motion indices
    i_gm = np.argwhere( df_prc_info.station == n_s ).flatten()
    #re-order components
    i_gm = i_gm[[np.argwhere(df_prc_info.loc[i_gm,'cmp'] == c)[0][0] for c in cmp_names]]
    assert( np.all(df_prc_info.loc[i_gm,'cmp'] == cmp_names) ),'Error. Incorrect order of component'
    
    #load time histories
    gm_acc = [pd.read_csv(dir_gm+df_prc_info.loc[i_gm[j],'fname_acc'], names=('time','processed','vel_pulse'), sep="\s+") 
              for j, c in enumerate(cmp_names)]
    gm_vel = [pd.read_csv(dir_gm+df_prc_info.loc[i_gm[j],'fname_vel'], names=('time','processed','vel_pulse'), sep="\s+") 
              for j, c in enumerate(cmp_names)]
    gm_dis = [pd.read_csv(dir_gm+df_prc_info.loc[i_gm[j],'fname_disp'], names=('time','processed','vel_pulse'), sep="\s+") 
              for j, c in enumerate(cmp_names)]

    #plot acceleration time histories
    # - - - - - - - - - - - - 
    for j1, (c, g_acc, g_vel, g_dis) in enumerate(zip(cmp_names,gm_acc,gm_vel,gm_dis)):
        #create figure
        # fig, axs = plt.subplots(3, figsize = (12,10))
        if j1 == 0: fig, axs = plt.subplots(3, figsize = (12,10))
        else: fig, axs = plt.subplots(3, figsize = (11.5,10))
        #plot acceleration time histories
        for j2, th in enumerate(['processed','vel_pulse']):
            lb = label[th]
            axs[0].plot(g_acc.time, g_acc[th], label=lb, 
                         color=lc_map[j2], linestyle=ls_map[j2], linewidth=lw_map[j2]) 
    
        #fig properties
        axs[0].grid()
        #labels
        axs[0].set_ylim(ylim_acc[k]*np.array([-1,1]))
        axs[0].yaxis.set_major_formatter('{x:.2f}')
        axs[0].tick_params(labelsize=25) 
        # axs[0].set_ylabel(c+f' - Acc.\n($m/sec^2$)',   fontsize=30)
        if j1==0: axs[0].set_ylabel(f'Acc. ($m/sec^2$)',   fontsize=30)
        #legend
        if j1 == 0: axs[0].legend(loc='lower right', ncol=2, fontsize=25)
    
        #plot velocity time histories
        for j2, th in enumerate(['processed','vel_pulse']):
            lb = label[th]
            axs[1].plot(g_vel.time, g_vel[th], label=lb, 
                         color=lc_map[j2], linestyle=ls_map[j2], linewidth=lw_map[j2]) 

        #fig properties
        axs[1].grid()
        #labels
        axs[1].set_ylim(ylim_vel[k]*np.array([-1,1]))
        axs[1].yaxis.set_major_formatter('{x:.0f}')
        axs[1].tick_params(labelsize=25) 
        # axs[1].set_ylabel(c+f' - Vel.\n($cm/sec$)',   fontsize=30)
        if j1==0: axs[1].set_ylabel(f'Vel. ($cm/sec$)',   fontsize=30)
        #plot displacement time histories
        for j2, th in enumerate(['processed','vel_pulse']):
            axs[2].plot(g_dis.time, g_dis[th], label=th, color=lc_map[j2], linestyle=ls_map[j2], linewidth=lw_map[j2]) 
    
        #fig properties
        axs[2].grid()
        #labels
        axs[2].set_ylim(ylim_dis[k]*np.array([-1,1]))
        axs[2].yaxis.set_major_formatter('{x:.0f}')
        axs[2].tick_params(labelsize=25) 
        # axs[2].set_ylabel(c+f' - Disp.\n($cm$)',   fontsize=30)
        if j1==0: axs[2].set_ylabel(f'Disp. ($cm$)',   fontsize=30)
        axs[2].set_xlabel(f'Time ($sec$)',   fontsize=30)
        #save time histories
        fig.tight_layout()
        fig.savefig( dir_fig+n_s+'_'+c+'_vel_pulse'+'.png', bbox_inches='tight')
        
    

