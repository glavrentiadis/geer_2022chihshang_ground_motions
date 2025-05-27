%% Summarize GM data for QGIS
%
%   Summarize PGA, PGV, P2PV, and permanent offset of processed ground motions
%   for gis comparison

%% Define Input
%flag input
flag_io = 2;
%flag roated
flag_rot = 0;

%tectonic offset
t_win = 1;
%peak to peak vel search window
p2p_win = 1.;

%components
if     flag_rot == 0; cmp = {'Z','N','E'};
elseif flag_rot == 1; cmp = {'Z','FN','FP'};
end
%horizontal compoents
i_horiz = [2,3];

%input files and directories
dir_inpt = '../../Data/ground_motions/';
if     flag_rot == 0; dir_inpt = [dir_inpt,'corrected_gm/'];
elseif flag_rot == 1; dir_inpt = [dir_inpt,'corrected_gm_rot/'];
end
if flag_io == 1
    fn_eq    = '2022_Guanshan';
    dir_inpt = [dir_inpt,'M6.5_0917/'];
elseif flag_io == 2 
    fn_eq    = '2022_Chihshang';
    dir_inpt = [dir_inpt,'M6.9_0918/'];
end

%output directories
dir_out =  '../../Data/gis/ground_motions/';
dir_fig = [dir_out,'figures/'];

%table options
warning('off','MATLAB:table:RowsAddedNewVars');

%% Load Data
%load ground motions
fn_df_gm = [fn_eq,'_gm_info_tilt_corrected.mat'];
load([dir_inpt,fn_df_gm ],'gm_prc_all','df_gm_info','df_gm_prc')
n_gm = size(df_gm_info,1);

%% Processing
%initialize 
df_gm_tect_info = df_gm_info;

%interate over ground motions
for k = 1:n_gm
    fprintf('Processing ground motion %i of %i ...\n',k,n_gm)
    %extract ground-motion records
    % ---   ---   ---   ---   ---   ---   ---
    %subset ground-motion structure
    gm = gm_prc_all(k,:);
    %subset ground motion components
    i_gm = all([strcmp(df_gm_prc.event,df_gm_tect_info{k,'event'}),     ...
                strcmp(df_gm_prc.network,df_gm_tect_info{k,'network'}), ...
                strcmp(df_gm_prc.station,df_gm_tect_info{k,'station'})], 2);
    df_gm = df_gm_prc(i_gm,:);

    %summarize ground-motion parameters
    % ---   ---   ---   ---   ---   ---   ---
    %time
    time = gm{1}.time;
    t_e = max(time);
    %acceleration
    acc = nan(length(time),length(cmp));
    for j = 1:length(gm)
        acc(:,strcmp(cmp,gm{j}.component)) = gm{j}.acc;
    end
    %velocity 
    vel = nan(length(time),length(cmp));
    for j = 1:length(gm)
        vel(:,strcmp(cmp,gm{j}.component)) = gm{j}.vel;
    end
    %displacement
    dis = nan(length(time),length(cmp));
    for j = 1:length(gm)
        dis(:,strcmp(cmp,gm{j}.component)) = gm{j}.dis;
    end

    %sumarize PGA
    % - - - - - - - - - - - - - - - -
    f_pga = cellfun(@(c) ['PGA_',c], cmp,'UniformOutput',false);
    [pga,i_pga] = max(abs(acc),[],1);
    df_gm_tect_info{k,f_pga} = diag(acc(i_pga,:))';
    %total horizontal component
    df_gm_tect_info{k,'PGA_horiz'} = norm(pga(i_horiz));
    %sumarize PGAi (instatenious time)
    % - - - - - - - - - - - - - - - -
    f_pgai = cellfun(@(c) ['PGAi_',c], cmp(i_horiz),'UniformOutput',false);
    [pgai,i_pgai] = max(vecnorm(acc(:,i_horiz) ,2,2)); %find pgai instance
    df_gm_tect_info{k,f_pgai} = acc(i_pgai,i_horiz);
    %total horizontal component
    df_gm_tect_info{k,'PGAi_horiz'} = pgai;

    %sumarize PGV
    % - - - - - - - - - - - - - - - -
    f_pgv = cellfun(@(c) ['PGV_',c], cmp,'UniformOutput',false);
    [pgv,i_pgv] = max(abs(vel),[],1);
    df_gm_tect_info{k,f_pgv} = diag(vel(i_pgv,:))';
    %total horizontal component
    df_gm_tect_info{k,'PGV_horiz'} = norm(pgv(i_horiz));
    %sumarize PGVi (instatenious time)
    % - - - - - - - - - - - - - - - -
    f_pgvi = cellfun(@(c) ['PGVi_',c], cmp(i_horiz),'UniformOutput',false);
    [pgvi,i_pgvi] = max(vecnorm(vel(:,i_horiz) ,2,2)); %find pgai instance
    df_gm_tect_info{k,f_pgvi} = vel(i_pgvi,i_horiz);
    %total horizontal component
    df_gm_tect_info{k,'PGVi_horiz'} = pgvi;

    %rolling peak to peak velocity 
    % - - - - - - - - - - - - - - - -
    f_ppv = cellfun(@(c) ['PPVi_',c], cmp(i_horiz),'UniformOutput',false);
    %initalize
    ppv = 0;
    df_gm_tect_info{k,f_ppv} = 0;
    %iterate widening peak to peak window up to p2p_win sec
    for j = reshape(find(time<=p2p_win),1,[]) %iterate over all time delays
        p2pv_horiz = vel(j:end,i_horiz) - vel(1:end-j+1,i_horiz);
        if max(vecnorm(p2pv_horiz,2,2)) > ppv
            [ppv,i_ppv] = max(vecnorm(p2pv_horiz,2,2)); 
            df_gm_tect_info{k,f_ppv} = p2pv_horiz(i_ppv,:);
            %total horizontal component
            df_gm_tect_info{k,'PPVi_horiz'} = ppv;
            %pulse period
            df_gm_tect_info{k,'PPVi_T'}  = 2 * time(j);
            df_gm_tect_info{k,'PPVi_T1'} = time(i_ppv);
            df_gm_tect_info{k,'PPVi_T2'} = time(i_ppv+j);
        end
    end
    
    %tectonic offset
    % - - - - - - - - - - - - - - - -
    f_tslip = cellfun(@(c) ['tect_slip_',c], cmp,'UniformOutput',false);
    %final displacement
    d_final = mean(dis(time >= t_e-t_win,:),1);
    df_gm_tect_info{k,f_tslip} = d_final;
    %total horizontal component
    df_gm_tect_info{k,'tect_slip_horiz'} = norm(df_gm_tect_info{k,f_tslip(i_horiz)});
end


%% Output
%create directories
if not(isfolder(dir_out)); mkdir(dir_out); end

%save tectonic offset
if     flag_rot == 0; fn_summary = [fn_eq,'_gm_info_summary.csv'];
elseif flag_rot == 1; fn_summary = [fn_eq,'_gm_info_summary_rot.csv'];
end
writetable(df_gm_tect_info,[dir_out,fn_summary])
