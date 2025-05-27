#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 21 19:04:32 2025

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
#input filenames
fn_vel_pulse = ['2022_Chihshangvel_pulse_info.csv','2022_Guanshanvel_pulse_info.csv']

#data directory
dir_data = '../../../Data/ground_motions/vel_pulse/summary/'

#output directory
dir_out = dir_data
dir_fig = dir_out + 'figures/'

#read velocity pulse data
df_vel_pulse = pd.concat([pd.read_csv(dir_data + fn_v_p) for fn_v_p in fn_vel_pulse], ignore_index=True) 

#
i_H  = np.isin(df_vel_pulse.loc[:,'cmp'],['N','E'])
i_V  = np.isin(df_vel_pulse.loc[:,'cmp'],['Z','V'])
#fling and directivity records
i_Ff = df_vel_pulse.loc[:,'F_fling'] == 1
i_Fd = df_vel_pulse.loc[:,'F_direct'] == 1


# %% Summary
# ======================================
msg = list()
msg.append(f'Summary compiled stong-motions:')
msg.append(f'\tFling Pulse Records: %i'%(sum(i_Ff)))
msg.append(f'\t\tHorizontal Component: %i'%(sum(np.logical_and(i_Ff,i_H))))
msg.append(f'\t\tVertical Component: %i'%(sum(np.logical_and(i_Ff,i_V))))
msg.append(f'\tDynamic Pulse Records: %i'%(sum(i_Fd)))
msg.append(f'\t\tHorizontal Component: %i'%(sum(np.logical_and(i_Fd,i_H))))
msg.append(f'\t\tVertical Component: %i'%(sum(np.logical_and(i_Fd,i_V))))

#print summary
for m in msg:
    print(m)

# %% Plotting
# ======================================
#create output directories
if not os.path.isdir(dir_fig): pathlib.Path(dir_fig).mkdir(parents=True, exist_ok=True)

#clear pulse timing histogram
fig, ax = plt.subplots(figsize=(10, 10))
fname_fig = 'pulse_timing'
#histogram for fling step
ax.hist(df_vel_pulse.loc[i_Ff, 'IvNorm_t1'], bins=np.arange(0, 1.1, 0.1), 
        weights=np.ones(sum(i_Ff)) / sum(i_Ff),  alpha=0.70,
        label='Fling-Step Pulse')
#histogram of directivity
ax.hist(df_vel_pulse.loc[i_Fd, 'IvNorm_t1'], bins=np.arange(0, 1.1, 0.1), 
        weights=np.ones(sum(i_Fd)) / sum(i_Fd),  alpha=0.7,
        label='Dynamic Pulse')
#figure properties
ax.set_xlabel("Normalized cumulative velocity, $i_v(t'_p)$", fontsize=30)
ax.set_ylabel('Probability',               fontsize=30)
ax.grid(which='both')
ax.set_ylim([0, 1])
ax.tick_params(axis='x', labelsize=28)
ax.tick_params(axis='y', labelsize=28)
ax.legend(fontsize=30)
fig.savefig( dir_fig + fname_fig + '.png', dpi=300, bbox_inches='tight') 

#clear pulse timing histogram
fig, ax = plt.subplots(figsize=(10, 10))
fname_fig = 'pulse_frequency'
#velocity pulse period
ax.hist(df_vel_pulse.loc[i_Fd, 'fp'], color='C1',
        weights=np.ones(sum(i_Fd)) / sum(i_Fd), bins=np.arange(0, 1.1, 0.1))
#figure properties
ax.set_xlabel('Pulse frequency, $f_p$ (hz)', fontsize=30)
ax.set_ylabel('Probability',             fontsize=30)
ax.grid(which='both')
ax.set_ylim([0, 1])
ax.tick_params(axis='x', labelsize=28)
ax.tick_params(axis='y', labelsize=28)
# ax.legend(fontsize=30)
fig.savefig( dir_fig + fname_fig + '.png', dpi=300, bbox_inches='tight') 
