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
#flag bs method
flag_bs = True

#gravity
gravity = 9.80665

#components
cmp_names = ['N','E','Z']

#station information
n_sta  = 'TSMIP - HWA075'
fn_sta = 'TSMIP_HWA075'

#ground motion data
dir_gm = '../../../Data/ground_motions/corrected_gm/example/' 
dir_gm += fn_sta
dir_gm += '_M6.9_0918_' + ('correct' if flag_bs else 'wrong') + '/'

#output directory
dir_out = '../../../Data/ground_motions/summary_general/'
dir_fig = dir_out + 'figures/'

# %% Plotting
# ======================================
#create output directories
if not os.path.isdir(dir_out): pathlib.Path(dir_out).mkdir(parents=True, exist_ok=True)
if not os.path.isdir(dir_fig): pathlib.Path(dir_fig).mkdir(parents=True, exist_ok=True)

#list gm file names
fn_gm_all = os.listdir(dir_gm)

#station filename
i_s = np.argwhere( np.array([bool(re.match('.*%s_%s\.acc'%(fn_sta,cmp_names[0]), f_gm)) for f_gm in fn_gm_all]) )[0][0]
f_s = re.search('(.*%s)_%s\.acc'%(fn_sta,cmp_names[0]), fn_gm_all[i_s]).group(1) 

#load time histories
gm_acc = [pd.read_csv('%s%s_%s.acc'%(dir_gm,f_s,c), names=('time','data'),  sep="\s+") for c in cmp_names]
gm_vel = [pd.read_csv('%s%s_%s.vel'%(dir_gm,f_s,c), names=('time','data'),  sep="\s+") for c in cmp_names]
gm_dis = [pd.read_csv('%s%s_%s.disp'%(dir_gm,f_s,c), names=('time','data'), sep="\s+") for c in cmp_names]

    

for k, (f_s, n_s) in enumerate(zip(fn_sta, n_sta)):
    #station filename
    i_s = np.argwhere( np.array([bool(re.match('.*%s_%s\.acc'%(f_s,cmp_names[0]), f_gm)) for f_gm in fn_gm_all]) )[0][0]
    f_s = re.search('(.*%s)_%s\.acc'%(f_s,cmp_names[0]), fn_gm_all[i_s]).group(1) 
   
    #load time histories
    gm_acc = [pd.read_csv('%s%s_%s.acc'%(dir_gm,f_s,c), names=('time','data'),  sep="\s+") for c in cmp_names]
    gm_vel = [pd.read_csv('%s%s_%s.vel'%(dir_gm,f_s,c), names=('time','data'),  sep="\s+") for c in cmp_names]
    gm_dis = [pd.read_csv('%s%s_%s.disp'%(dir_gm,f_s,c), names=('time','data'), sep="\s+") for c in cmp_names]
    
    #create figure
    fig, axs = plt.subplots(3, figsize = (12,10))
    #plot acceleration time histories
    for g_acc, cmp in zip(gm_acc, cmp_names):
        axs[0].plot(g_acc.time, g_acc.data/gravity, linewidth=2, label=cmp)    
    #fig properties
    axs[0].grid()
    if k==0: axs[0].legend(loc='upper right', fontsize=35)
    #labels
    axs[0].set_ylim([-0.6,0.6])
    axs[1].yaxis.set_major_formatter('{x:.2f}')
    axs[0].tick_params(labelsize=25) 
    axs[0].set_xticklabels([])
    axs[0].set_ylabel('Acc. (g)',   fontsize=30)
    # axs[0].set_xlabel('Time (sec)', fontsize=30)
    #plot velocity time histories
    for g_vel, cmp in zip(gm_vel, cmp_names):
        axs[1].plot(g_vel.time, g_vel.data, linewidth=2, label=cmp)    
    #fig properties
    axs[1].grid()
    #labels
    axs[1].set_ylim([-110,110])
    axs[1].yaxis.set_major_formatter('{x:.0f}')
    axs[1].tick_params(labelsize=25) 
    axs[1].set_ylabel('Vel. (cm/sec)',   fontsize=30)
    # axs[1].set_xlabel('Time (sec)', fontsize=30)
    #plot velocity time histories
    for g_dis, cmp in zip(gm_dis, cmp_names):
        axs[2].plot(g_dis.time, g_dis.data, linewidth=2, label=cmp)    
    #fig properties
    axs[2].grid()
    #labels
    axs[2].set_ylim([-110,110])
    axs[2].yaxis.set_major_formatter('{x:.0f}')
    axs[2].tick_params(labelsize=25) 
    axs[2].set_ylabel('Disp. (cm)',   fontsize=30)
    axs[2].set_xlabel('Time (sec)', fontsize=30)
    
    
    #save time histories
    fig.tight_layout()
    fig.savefig( dir_fig+f_s+'.png', bbox_inches='tight')
    
