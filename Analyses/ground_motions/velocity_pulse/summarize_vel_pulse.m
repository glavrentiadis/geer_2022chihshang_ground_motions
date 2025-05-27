%% Sumarize Velocity Pulses

%

%% Define Input
addpath('matlab_lib/')
addpath('matlab_lib/ground_motions/')
addpath('matlab_lib/plotting/')

%flag input
% flag_io=1: M6.5 2022 Guanshan earthquake
% flag_io=2: M6.9 2022 Chihshang earthquake
flag_io = 1;

%flag roated
flag_rot = 0;

%filtering parameters
%forward directivity
fd_dres_ratio = 0.5;
fd_vp_ratio   = 0.5;
fd_gamma      = 2.0;
%fling effect
fl_dres_ratio = 0.5;
fl_gamma      = 2.0;

% directories
% ---   ---   ---   ---   ---
%input files and directories
dir_inpt = '../../../Data/ground_motions/';
if     flag_rot == 0; dir_inpt = [dir_inpt,'corrected_gm/'];
elseif flag_rot == 1; dir_inpt = [dir_inpt,'corrected_gm_rot/'];
end

%earthquake file
if     flag_io == 1
    fn_eq    = '2022_Guanshan';
    dir_inpt = [dir_inpt,'M6.5_0917/'];
elseif flag_io == 2
    fn_eq    = '2022_Chihshang';
    dir_inpt = [dir_inpt,'M6.9_0918/'];
end

dir_data = '../../../Data/ground_motions/';
if     flag_rot == 0; dir_data = [dir_data,'vel_pulse/'];
elseif flag_rot == 1; dir_data = [dir_data,'vel_pulse_rot/'];
end
if     flag_io == 1; sdir_data = 'M6.5_0917/';
elseif flag_io == 2; sdir_data = 'M6.9_0918/';
end

%ground-motion component information
fn_gm_info = [fn_eq,'_gm_info_vel_pulse.mat'];

%output directories
dir_out = dir_data;
dir_fig = [dir_out,'figures/'];

%% Read Input
%load ground motion data
load([dir_data,sdir_data,fn_gm_info],'df_gm_info','df_vel_pulse','gm_inst_all','vel_pulse_all')

%% Processing
%velocity pulse ratio
df_vel_pulse{:,'vp_dres_ratio'} = df_vel_pulse{:,'vp_dres'} ./ df_vel_pulse{:,'vp_dpeak'};

i = 1;
for k = 1:size(gm_inst_all,1)
    for j = 1:size(gm_inst_all,2)
        %velocity parameters
        t0        = df_vel_pulse{i,{'t0'}};
        A         = df_vel_pulse{i,{'A'}};
        fp        = df_vel_pulse{i,{'fp'}};
        nu        = df_vel_pulse{i,{'nu'}};
        gamma     = df_vel_pulse{i,{'gamma'}};
        %start time
        ts        = t0-gamma/(2*fp);
        t1        = t0-gamma/(4*fp);

        %velocity parameters
        pgv  = max(abs(gm_inst_all{k,j}.vel));
        vp   = max(abs(vel_pulse_all{k,j}.vel));
        %displacement parameters
        pgd  = max(abs(gm_inst_all{k,j}.dis));
        % dres = df_vel_pulse{i,'vp_dres'};
        dres = abs(vel_pulse_all{k,j}.dis(end));
        
        %identify pulse time
        [~,idx_ts] = min(abs(gm_inst_all{k,j}.time - ts));
        [~,idx_t1] = min(abs(gm_inst_all{k,j}.time - t1));

        %normalized arias velocity at velocity pulse
        df_vel_pulse{i,{'IvNorm_ts'}} = gm_inst_all{k,j}.IvNorm_his(idx_ts);
        df_vel_pulse{i,{'IvNorm_t0'}} = gm_inst_all{k,j}.IvNorm_his(idx_t0);
        df_vel_pulse{i,{'IvNorm_t1'}} = gm_inst_all{k,j}.IvNorm_his(idx_t1);
        
        %summarize velocity parameters
        df_vel_pulse{i,{'pgv'}}        = pgv;
        df_vel_pulse{i,{'vp'}}         = vp;
        df_vel_pulse{i,{'vp_ratio'}}   = vp / pgv;
        %summarize displacement parameters 
        df_vel_pulse{i,{'pgd'}}        = pgd;
        df_vel_pulse{i,{'dres'}}       = dres;
        df_vel_pulse{i,{'dres_ratio'}} = dres / pgd;        

        %index 
        i = i + 1;
    end
end

%identify directivity and fling and velocity pulses
i_fd = df_vel_pulse.gamma < fd_gamma & df_vel_pulse.vp_dres_ratio < fd_dres_ratio & df_vel_pulse.vp_ratio > fd_vp_ratio;
i_fl = df_vel_pulse.gamma < fl_gamma & df_vel_pulse.dres_ratio > fl_dres_ratio;
%store record classification
df_vel_pulse{:,'F_direct'} = i_fd;
df_vel_pulse{:,'F_fling'} = i_fl;

%% Output
%create directories
if not(isfolder(dir_out)); mkdir(dir_out); end
if not(isfolder(dir_fig)); mkdir(dir_fig); end

%writeout vel pulse info
writetable(df_vel_pulse,[dir_out,fn_eq,'vel_pulse_info','.csv'])

