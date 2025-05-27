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
flag_bs_corr = True

#flag annotation
flag_annot = False

#gravity
gravity = 9.80665

#station information
n_sta  = 'TSMIP - HWA075'
fn_sta = 'HWA075'

#components
cmp_names = ['N','E','Z']

#ground motion data
dir_gm = '../../../Data/ground_motions/example_timing/' 
dir_gm += 'TSMIP_' + fn_sta + '_M6.9_0918/'

#output directory
dir_out = dir_gm
dir_fig = dir_out + 'figures/'

#axis limits
ylim_acc = 6.0
ylim_vel = 75.
ylim_dis = 50.

#line parameters
lw_map = [2,2,1]
lc_map = ['#ff7f0e','#1f77b4','black']
ls_map  = ['-','-','--']

#ground motion labels
label = {'seed':'Original',
         'baseline corrected':'Baseline Corrected',
         'correction':None}

#filename processing info
fn_prc_info  = '2022_Chihshang_gm_info_tilt_corrected_examp_'
fn_prc_info += ('correct' if flag_bs_corr else 'wrong') + '.csv'

# %% Plotting
# ======================================
#create output directories
if not os.path.isdir(dir_fig): pathlib.Path(dir_fig).mkdir(parents=True, exist_ok=True)

#output filename
fn_sta_out = fn_sta + ('_correct' if flag_bs_corr else '_wrong')

#processing info
df_prc_info = pd.read_csv(dir_gm+fn_prc_info)

#read ground motions
# - - - - - - - - - - - - 
#ground motion indices
i_gm = np.argwhere( df_prc_info.station == fn_sta ).flatten()
#re-order components
i_gm = i_gm[[np.argwhere(df_prc_info.loc[i_gm,'component'] == c)[0][0] for c in cmp_names]]
assert( np.all(df_prc_info.loc[i_gm,'component'] == cmp_names) ),'Error. Incorrect order of component'

#load time histories
gm_acc = [pd.read_csv(dir_gm+df_prc_info.loc[i_gm[j],'fname_acc'], names=('time','baseline corrected','seed','correction'), sep="\s+") 
          for j, c in enumerate(cmp_names)]
gm_vel = [pd.read_csv(dir_gm+df_prc_info.loc[i_gm[j],'fname_vel'], names=('time','baseline corrected','seed','correction'), sep="\s+") 
          for j, c in enumerate(cmp_names)]
gm_dis = [pd.read_csv(dir_gm+df_prc_info.loc[i_gm[j],'fname_disp'], names=('time','baseline corrected','seed','correction'), sep="\s+") 
          for j, c in enumerate(cmp_names)]

#plot acceleration time histories
# - - - - - - - - - - - - 
#create figure
fig, axs = plt.subplots(3, figsize = (12,10))
for j1, (c,g_acc) in enumerate( zip(cmp_names,gm_acc) ) :
    #plot velocity time histories
    for j2, th in enumerate(['seed','baseline corrected','correction']):
        lb = label[th]
        axs[j1].plot(g_acc.time, g_acc[th], label=lb, 
                     color=lc_map[j2], linestyle=ls_map[j2], linewidth=lw_map[j2]) 

    #base line correction info
    prc_info = df_prc_info.loc[i_gm[j1],:]

    if prc_info.method == 'BoorePieceWise':
        axs[j1].axvline(x=prc_info.t1, color='gray', linestyle='--', linewidth=1)
        axs[j1].axvline(x=prc_info.t2, color='gray', linestyle='--', linewidth=1)
        if j1==0:
            axs[j1].annotate('$t_1$', xy = (prc_info.t1, 0), xytext = (prc_info.t1-5, .5),
                             arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                             horizontalalignment = 'center',  fontsize=25)
            
            i_t2 = np.argmin(np.abs(g_acc.time-prc_info.t2))
            a_t2 =  g_acc.loc[i_t2,'correction']
            axs[j1].annotate('$t_2$', xy = (prc_info.t2, a_t2), xytext = (prc_info.t2+5, (abs(a_t2)+.5)),
                             arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                             horizontalalignment = 'center',  fontsize=25)            
    elif prc_info.method == 'BooreQuad':
            axs[j1].axvline(x=prc_info.t1, color='gray', linestyle='--', linewidth=1)
            if j1==0:
                axs[j1].annotate('$t_q$', xy = (prc_info.t1, 0), xytext = (prc_info.t1-5, .5),
                                 arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                                 horizontalalignment = 'center',  fontsize=25)
    #fig properties
    axs[j1].grid()
    #labels
    axs[j1].set_ylim(ylim_acc*np.array([-1,1]))
    axs[j1].yaxis.set_major_formatter('{x:.2f}')
    axs[j1].tick_params(labelsize=25) 
    axs[j1].set_ylabel(c+f' - Acc.\n($m/sec^2$)',   fontsize=30)
    #legend
    if j1 == 0: axs[j1].legend(loc='upper right', fontsize=25)
    
#save time histories
fig.tight_layout()
fig.savefig( dir_fig+fn_sta_out+'_acc'+'.png', bbox_inches='tight')
    
#plot acceleration time histories (without seed and correction)
# - - - - - - - - - - - - 
#create figure
fig, axs = plt.subplots(3, figsize = (12,10))
for j1, (c,g_acc) in enumerate( zip(cmp_names,gm_acc) ) :
    #plot velocity time histories
    j2, th = 1, 'baseline corrected'
    axs[j1].plot(g_acc.time, g_acc[th], label=th, color=lc_map[j2], linestyle=ls_map[j2], linewidth=lw_map[j2]) 

    #base line correction info
    prc_info = df_prc_info.loc[i_gm[j1],:]
    
    if prc_info.method == 'BoorePieceWise':
        axs[j1].axvline(x=prc_info.t1, color='gray', linestyle='--', linewidth=1)
        axs[j1].axvline(x=prc_info.t2, color='gray', linestyle='--', linewidth=1)
        if j1==0:
            axs[j1].annotate('$t_1$', xy = (prc_info.t1, 0), xytext = (prc_info.t1-5, .5),
                             arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                             horizontalalignment = 'center',  fontsize=25)
            
            i_t2 = np.argmin(np.abs(g_acc.time-prc_info.t2))
            a_t2 =  0. * g_acc.loc[i_t2,'correction']
            axs[j1].annotate('$t_2$', xy = (prc_info.t2, a_t2), xytext = (prc_info.t2+5, (abs(a_t2)+.5)),
                             arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                             horizontalalignment = 'center',  fontsize=25)            
    elif prc_info.method == 'BooreQuad':
            axs[j1].axvline(x=prc_info.t1, color='gray', linestyle='--', linewidth=1)
            if j1==0:
                axs[j1].annotate('$t_q$', xy = (prc_info.t1, 0), xytext = (prc_info.t1-5, .5),
                                 arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                                 horizontalalignment = 'center',  fontsize=25)
    #fig properties
    axs[j1].grid()
    #labels            
    axs[j1].set_ylim(ylim_acc*np.array([-1,1]))
    axs[j1].yaxis.set_major_formatter('{x:.2f}')
    axs[j1].tick_params(labelsize=25) 
    axs[j1].set_ylabel(c+f' - Acc.\n($m/sec^2$)',   fontsize=30)
    
#save time histories
fig.tight_layout()
fig.savefig( dir_fig+fn_sta_out+'_acc_corr'+'.png', bbox_inches='tight')

#plot velocity time histories
# - - - - - - - - - - - - 
#create figure
fig, axs = plt.subplots(3, figsize = (12,10))
for j1, (c,g_vel) in enumerate( zip(cmp_names,gm_vel) ) :
    #plot velocity time histories
    for j2, th in enumerate(['seed','baseline corrected','correction']):
        lb = label[th]
        axs[j1].plot(g_vel.time, g_vel[th], label=lb, 
                     color=lc_map[j2], linestyle=ls_map[j2], linewidth=lw_map[j2]) 

    #base line correction info
    prc_info = df_prc_info.loc[i_gm[j1],:]

    if flag_annot:
        if prc_info.method == 'BoorePieceWise':
            axs[j1].axvline(x=prc_info.t1, color='gray', linestyle='--', linewidth=1)
            axs[j1].axvline(x=prc_info.t2, color='gray', linestyle='--', linewidth=1)
            if j1==0:
                axs[j1].annotate('$t_1$', xy = (prc_info.t1, 0), xytext = (prc_info.t1-5, 10),
                                 arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                                 horizontalalignment = 'center',  fontsize=25)
                
                i_t2 = np.argmin(np.abs(g_vel.time-prc_info.t2))
                v_t2 =  g_vel.loc[i_t2,'correction']
                axs[j1].annotate('$t_2$', xy = (prc_info.t2, v_t2), xytext = (prc_info.t2+5, (abs(v_t2)+10)),
                                 arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                                 horizontalalignment = 'center',  fontsize=25)            
        elif prc_info.method == 'BooreQuad':
                axs[j1].axvline(x=prc_info.t1, color='gray', linestyle='--', linewidth=1)
                if j1==0:
                    axs[j1].annotate('$t_q$', xy = (prc_info.t1, 0), xytext = (prc_info.t1-5, 10),
                                     arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                                     horizontalalignment = 'center',  fontsize=25)
    #fig properties
    axs[j1].grid()
    #labels
    axs[j1].set_ylim(ylim_vel*np.array([-1,1]))
    axs[j1].yaxis.set_major_formatter('{x:.0f}')
    axs[j1].tick_params(labelsize=25) 
    axs[j1].set_ylabel(c+f' - Vel.\n($cm/sec$)',   fontsize=30)
    #legend
    if j1 == 0: axs[j1].legend(loc='lower right', fontsize=25)
    
#save time histories
fig.tight_layout()
fig.savefig( dir_fig+fn_sta_out+'_vel'+'.png', bbox_inches='tight')

#plot velocity time histories (without seed and correction)
# - - - - - - - - - - - - 
#create figure
fig, axs = plt.subplots(3, figsize = (12,10))
for j1, (c,g_vel) in enumerate( zip(cmp_names,gm_vel) ) :
    #plot velocity time histories
    j2, th = 1, 'baseline corrected'
    lb = label[th]
    axs[j1].plot(g_vel.time, g_vel[th], label=lb, 
                 color=lc_map[j2], linestyle=ls_map[j2], linewidth=lw_map[j2]) 

    #base line correction info
    prc_info = df_prc_info.loc[i_gm[j1],:]

    if flag_annot:
        if prc_info.method == 'BoorePieceWise':
            axs[j1].axvline(x=prc_info.t1, color='gray', linestyle='--', linewidth=1)
            axs[j1].axvline(x=prc_info.t2, color='gray', linestyle='--', linewidth=1)
            if j1==0:
                axs[j1].annotate('$t_1$', xy = (prc_info.t1, 0), xytext = (prc_info.t1-5, 10),
                                 arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                                 horizontalalignment = 'center',  fontsize=25)
                
                i_t2 = np.argmin(np.abs(g_vel.time-prc_info.t2))
                v_t2 =  g_vel.loc[i_t2,'correction']
                axs[j1].annotate('$t_2$', xy = (prc_info.t2, v_t2), xytext = (prc_info.t2+5, (abs(v_t2)+10)),
                                 arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                                 horizontalalignment = 'center',  fontsize=25)            
        elif prc_info.method == 'BooreQuad':
                axs[j1].axvline(x=prc_info.t1, color='gray', linestyle='--', linewidth=1)
                if j1==0:
                    axs[j1].annotate('$t_q$', xy = (prc_info.t1, 0), xytext = (prc_info.t1-5, 10),
                                     arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                                     horizontalalignment = 'center',  fontsize=25)
    #fig properties
    axs[j1].grid()
    #labels
    axs[j1].set_ylim(ylim_vel*np.array([-1,1]))
    axs[j1].yaxis.set_major_formatter('{x:.0f}')
    axs[j1].tick_params(labelsize=25) 
    axs[j1].set_ylabel(c+f' - Vel.\n($cm/sec$)',   fontsize=30)
    
#save time histories
fig.tight_layout()
fig.savefig( dir_fig+fn_sta_out+'_vel_corr'+'.png', bbox_inches='tight')

#plot displacement time histories
# - - - - - - - - - - - - 
#create figure
fig, axs = plt.subplots(3, figsize = (12,10))
for j1, (c,g_dis) in enumerate( zip(cmp_names,gm_dis) ) :
    #plot velocity time histories
    for j2, th in enumerate(['seed','baseline corrected','correction']):
        axs[j1].plot(g_dis.time, g_dis[th], label=th, color=lc_map[j2], linestyle=ls_map[j2], linewidth=lw_map[j2]) 

    #base line correction info
    prc_info = df_prc_info.loc[i_gm[j1],:]

    if flag_annot:
        if prc_info.method == 'BoorePieceWise':
            axs[j1].axvline(x=prc_info.t1, color='gray', linestyle='--', linewidth=1)
            axs[j1].axvline(x=prc_info.t2, color='gray', linestyle='--', linewidth=1)
            if j1==0:
                axs[j1].annotate('$t_1$', xy = (prc_info.t1, 0), xytext = (prc_info.t1-5, 10),
                                 arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                                 horizontalalignment = 'center',  fontsize=25)
                
                i_t2 = np.argmin(np.abs(g_vel.time-prc_info.t2))
                d_t2 =  g_dis.loc[i_t2,'correction']
                axs[j1].annotate('$t_2$', xy = (prc_info.t2, d_t2), xytext = (prc_info.t2+5, (abs(d_t2)+10)),
                                 arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                                 horizontalalignment = 'center',  fontsize=25)            
        elif prc_info.method == 'BooreQuad':
                axs[j1].axvline(x=prc_info.t1, color='gray', linestyle='--', linewidth=1)
                if j1==0:
                    axs[j1].annotate('$t_q$', xy = (prc_info.t1, 0), xytext = (prc_info.t1-5, 5),
                                     arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                                     horizontalalignment = 'center',  fontsize=25)
    #fig properties
    axs[j1].grid()
    #labels
    axs[j1].set_ylim(ylim_dis*np.array([-1,1]))
    axs[j1].yaxis.set_major_formatter('{x:.0f}')
    axs[j1].tick_params(labelsize=25) 
    axs[j1].set_ylabel(c+f' - Disp.\n($cm$)',   fontsize=30)
    
#save time histories
fig.tight_layout()
fig.savefig( dir_fig+fn_sta_out+'_dis'+'.png', bbox_inches='tight')

#plot displacement time histories (without seed and correction)
# - - - - - - - - - - - - 
#create figure
fig, axs = plt.subplots(3, figsize = (12,10))
for j1, (c,g_dis) in enumerate( zip(cmp_names,gm_dis) ) :
    #plot velocity time histories
    j2, th = 1, 'baseline corrected'
    axs[j1].plot(g_dis.time, g_dis[th], label=th, color=lc_map[j2], linestyle=ls_map[j2], linewidth=lw_map[j2]) 

    #base line correction info
    prc_info = df_prc_info.loc[i_gm[j1],:]

    if flag_annot:
        if prc_info.method == 'BoorePieceWise':
            axs[j1].axvline(x=prc_info.t1, color='gray', linestyle='--', linewidth=1)
            axs[j1].axvline(x=prc_info.t2, color='gray', linestyle='--', linewidth=1)
            if j1==0:
                axs[j1].annotate('$t_1$', xy = (prc_info.t1, 0), xytext = (prc_info.t1-5, 10),
                                 arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                                 horizontalalignment = 'center',  fontsize=25)
                
                i_t2 = np.argmin(np.abs(g_vel.time-prc_info.t2))
                d_t2 =  0. * g_dis.loc[i_t2,'correction']
                axs[j1].annotate('$t_2$', xy = (prc_info.t2, d_t2), xytext = (prc_info.t2+5, (abs(d_t2)+10)),
                                 arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                                 horizontalalignment = 'center',  fontsize=25)            
        elif prc_info.method == 'BooreQuad':
                axs[j1].axvline(x=prc_info.t1, color='gray', linestyle='--', linewidth=1)
                if j1==0:
                    axs[j1].annotate('$t_q$', xy = (prc_info.t1, 0), xytext = (prc_info.t1-5, 5),
                                     arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
                                     horizontalalignment = 'center',  fontsize=25)
    #fig properties
    axs[j1].grid()
    #labels
    axs[j1].set_ylim(ylim_dis*np.array([-1,1]))
    axs[j1].yaxis.set_major_formatter('{x:.0f}')
    axs[j1].tick_params(labelsize=25) 
    axs[j1].set_ylabel(c+f' - Disp.\n($cm$)',   fontsize=30)
    
#save time histories
fig.tight_layout()
fig.savefig( dir_fig+fn_sta_out+'_dis_corr'+'.png', bbox_inches='tight')
