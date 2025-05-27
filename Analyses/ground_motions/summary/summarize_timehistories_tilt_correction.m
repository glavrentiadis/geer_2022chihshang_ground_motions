%% Summarize Baseline Correction
%
% This script summarizes the baseline corrected time histories including
% the original, baseline corrected, and correction time histories

%% Define Input
%flag bs timing
% flag_io = 1;
flag_io = 2;
%flag roated
flag_rot = 0;

%gravity
gravity = 9.80665;

%station info
fn_sta = {'TTN001','TTN002','TTN014','TTN015','TTN020','TTN025','TTN026','TTN045','TTN061','HWA004','HWA037','HWA073'};
% fn_sta = {'2718'};

%ground motion data
dir_gm = '../../../Data/ground_motions/';
if     flag_rot == 0; dir_gm = [dir_gm,'corrected_gm/'];
elseif flag_rot == 1; dir_gm = [dir_gm,'corrected_gm_rot/'];
end
if flag_io == 1
    fn_eq  = '2022_Guanshan';
    dir_gm = [dir_gm,'M6.5_0917/'];
elseif flag_io == 2
    fn_eq  = '2022_Chihshang';
    dir_gm = [dir_gm,'M6.9_0918/'];
elseif flag_io == 3
    fn_eq    = '2023_Pazarcık_Turkey';
    dir_gm = [dir_gm,'M7.8_0206/'];
end

%output directories
dir_out = '../../../Data/ground_motions/summary/';
if     flag_rot == 0; dir_out = [dir_out,'corrected_gm_examp/'];
elseif flag_rot == 1; dir_out = [dir_out,'corrected_gm_rot_examp/'];
end
if flag_io == 1;     dir_out = [dir_out,'M6.5_0917/'];
elseif flag_io == 2; dir_out = [dir_out,'M6.9_0918/'];
elseif flag_io == 3; dir_out = [dir_out,'M7.8_0206/'];
end

%% Load Data
%load ground motions
fn_df_gm = [fn_eq,'_gm_info_tilt_corrected.mat'];
load([dir_gm,fn_df_gm ],'df_gm_info','df_gm_prc','gm_prc_all','gm_seed_all','gm_corr_all')

%% Processing
%create directories
if not(isfolder(dir_out)); mkdir(dir_out); end

%iterate over stations
n_sta = length(fn_sta);
for k = 1:n_sta
    fprintf('Processed ground motion: %s ( %i of %i)\n',fn_sta{k},k,n_sta)

    %select station
    i_gm  = find( strcmp(df_gm_prc.station,  fn_sta(k)) );
    i_sta = find( strcmp(df_gm_info.station, fn_sta(k)) ); 
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
        acc = array2table([time,prc_acc(:,j),seed_acc(:,j),corr_acc(:,j)], ...
                          'VariableNames',{'time','acc_prc','acc_seed','acc_corr'});
        writetable(acc,[dir_out,fname_acc],'Delimiter',' ','QuoteStrings','none','WriteVariableNames',false,'FileType','text')
        %velocity
        fname_vel = df_gm_prc{i_gm(j),'fname_vel'}{1};
        vel = array2table([time,prc_vel(:,j),seed_vel(:,j),corr_vel(:,j)], ...
                          'VariableNames',{'time','vel_prc','vel_seed','vel_corr'});
        writetable(vel,[dir_out,fname_vel],'Delimiter',' ','QuoteStrings','none','WriteVariableNames',false,'FileType','text')
        %displacement
        fname_dis = df_gm_prc{i_gm(j),'fname_disp'}{1};
        dis = array2table([time,prc_dis(:,j),seed_dis(:,j),corr_dis(:,j)], ...
                          'VariableNames',{'time','dis_prc','dis_seed','dis_corr'});
        writetable(dis,[dir_out,fname_dis],'Delimiter',' ','QuoteStrings','none','WriteVariableNames',false,'FileType','text')
    end
end

%save gm processing info
i_gm = reshape(cell2mat( cellfun(@(fn_s) find(strcmp(df_gm_prc.station,fn_s)), fn_sta,'UniformOutput',false) ),[],1);
fn_prc_info = [fn_eq,'_gm_info_tilt_corrected_examp.csv'];
writetable(df_gm_prc(i_gm,:),[dir_out,fn_prc_info])
