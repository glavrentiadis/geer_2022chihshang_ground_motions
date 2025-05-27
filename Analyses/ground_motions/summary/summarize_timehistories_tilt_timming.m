%% Summarize Baseline Correction (Timing Comparison)
%
% This script summarizes the baseline corrected time histories including
% the original, baseline corrected, and correction time histories

%% Define Input
%flag bs timing
flag_bs_corr = true;

%gravity
gravity = 9.80665;

%station info
fn_sta = 'HWA075';

%earthquake info
fn_eq  = '2022_Chihshang';

%ground motion data
dir_gm = '../../../Data/ground_motions/example_timing/';
dir_gm = [dir_gm, 'TSMIP_', fn_sta, '_M6.9_0918/'];
%subdirectory
if flag_bs_corr; sdir_gm = 'correct'; else; sdir_gm = 'wrong'; end

%output directories
dir_out = dir_gm;

%% Load Data
%load ground motions
fn_df_gm = [fn_eq,'_gm_info_tilt_corrected.mat'];
load([dir_gm,sdir_gm,'/',fn_df_gm ],'df_gm_info','df_gm_prc','gm_prc_all','gm_seed_all','gm_corr_all')

%% Processing
%create directories
if not(isfolder(dir_out)); mkdir(dir_out); end

%ground motion information
fprintf('Processed ground motion: %s \n',fn_sta)

%select station
i_gm  = find( strcmp(df_gm_prc.station,  fn_sta) );
i_sta = find( strcmp(df_gm_info.station, fn_sta) ); 
assert(~isempty(i_gm),'Error. Unavailable station.')

%ground motion info
gm_prc  = gm_prc_all(i_sta,:);
gm_seed = gm_seed_all(i_sta,:);
gm_corr = gm_corr_all(i_sta,:);
    
%summarize ground motion
time = gm_prc{1}.time;
prc_acc  = nan(length(time),3); prc_vel = nan(length(time),3);  prc_dis = nan(length(time),3);
seed_acc = nan(length(time),3); seed_vel = nan(length(time),3); seed_dis = nan(length(time),3);
corr_acc = nan(length(time),3); corr_vel = nan(length(time),3); corr_dis = nan(length(time),3);
for j = 1:length(gm_prc)
    prc_acc(:,j)  = gm_prc{j}.acc;  prc_vel(:,j)  = gm_prc{j}.vel;  prc_dis(:,j)  = gm_prc{j}.dis;
    seed_acc(:,j) = gm_seed{j}.acc; seed_vel(:,j) = gm_seed{j}.vel; seed_dis(:,j) = gm_seed{j}.dis;
    corr_acc(:,j) = gm_corr{j}.acc; corr_vel(:,j) = gm_corr{j}.vel; corr_dis(:,j) = gm_corr{j}.dis;
end
    
%save time histories
for j = 1:length(gm_prc)
    %acceleration
    fname_acc = df_gm_prc{i_gm(j),'fname_acc'}{1};
    fname_acc = replace(fname_acc,'.',['_',sdir_gm,'.']);
    acc = array2table([time,prc_acc(:,j),seed_acc(:,j),corr_acc(:,j)], ...
                      'VariableNames',{'time','acc_prc','acc_seed','acc_corr'});
    writetable(acc,[dir_out,fname_acc],'Delimiter',' ','QuoteStrings','none','WriteVariableNames',false,'FileType','text')
    %velocity
    fname_vel = df_gm_prc{i_gm(j),'fname_vel'}{1};
    fname_vel = replace(fname_vel,'.',['_',sdir_gm,'.']);
    vel = array2table([time,prc_vel(:,j),seed_vel(:,j),corr_vel(:,j)], ...
                      'VariableNames',{'time','vel_prc','vel_seed','vel_corr'});
    writetable(vel,[dir_out,fname_vel],'Delimiter',' ','QuoteStrings','none','WriteVariableNames',false,'FileType','text')
    %displacement
    fname_dis = df_gm_prc{i_gm(j),'fname_disp'}{1};
    fname_dis = replace(fname_dis,'.',['_',sdir_gm,'.']);
    dis = array2table([time,prc_dis(:,j),seed_dis(:,j),corr_dis(:,j)], ...
                      'VariableNames',{'time','dis_prc','dis_seed','dis_corr'});
    writetable(dis,[dir_out,fname_dis],'Delimiter',' ','QuoteStrings','none','WriteVariableNames',false,'FileType','text')

    %update filenames
    df_gm_prc{i_gm(j),'fname_acc'}{1}  = fname_acc;
    df_gm_prc{i_gm(j),'fname_vel'}{1}  = fname_vel;
    df_gm_prc{i_gm(j),'fname_disp'}{1} = fname_dis;
end

%save gm processing info
i_gm = find(strcmp(df_gm_prc.station,fn_sta));
fn_prc_info = [fn_eq,'_gm_info_tilt_corrected_examp_',sdir_gm,'.csv'];
writetable(df_gm_prc(i_gm,:),[dir_out,fn_prc_info])
