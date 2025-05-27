%% Process Accelearion Timehistories For Tilt Correction
%
%   This script perfoms baseline correction on raw acceleration time histories to remove the
%   the residual tilt from the ground motion

%% Define Input
addpath('matlab_lib/')
addpath('matlab_lib/ground_motions/')
addpath('matlab_lib/plotting/')
%flag input
flag_io = 2;
%flag roated
flag_rot = 0;

%acceleration uits
acc_grav = 9.80665;
acc_unt  = 'm/sec^2';

%components
if     flag_rot == 0; cmp = {'Z','N','E'};
elseif flag_rot == 1; cmp = {'Z','FN','FP'};
end

%psa periods
per = logspace(-2,1,300);

%input files and directories
dir_inpt = '../../../Data/ground_motions/';
if     flag_rot == 1; dir_inpt = [dir_inpt,'rotated_motions/'];
elseif flag_rot == 0; dir_inpt = [dir_inpt,'seed_motions/'];
end
if flag_io == 1
    fn_eq    = '2022_Guanshan';
    dir_inpt = [dir_inpt,'M6.5_0917/'];
elseif flag_io == 2 
    fn_eq    = '2022_Chihshang';
    dir_inpt = [dir_inpt,'M6.9_0918/'];
elseif flag_io == 3
    fn_eq    = '2023_Pazarcık_Turkey';
    dir_inpt = [dir_inpt,'M7.8_0206/'];
end

%output directories
dir_out = '../../../Data/ground_motions/';
if     flag_rot == 0; dir_out = [dir_out,'corrected_gm/'];
elseif flag_rot == 1; dir_out = [dir_out,'corrected_gm_rot/'];
end
if     flag_io == 1; dir_out = [dir_out,'M6.5_0917/'];
elseif flag_io == 2; dir_out = [dir_out,'M6.9_0918/'];
elseif flag_io == 3; dir_out = [dir_out,'M7.8_0206/'];
end
dir_fig = [dir_out,'figures/'];

%figure options
set (0,'DefaultFigurePaperType','usletter');
set (0,'DefaultFigureWindowStyle','docked');
%table options
warning('off','MATLAB:table:RowsAddedNewVars');


%% Read Input
%create directories
if not(isfolder(dir_out)); mkdir(dir_out); end
if not(isfolder(dir_fig)); mkdir(dir_fig); end

%load ground motions
fn_gm_prcd = [fn_eq,'_gm_info_tilt_corrected.mat'];
if exist([dir_out,fn_gm_prcd], 'file') == 2
    %load processed ground motions
    load([dir_out,fn_gm_prcd ],'df_gm_info','df_gm_prc','gm_prc_all','gm_seed_all','gm_corr_all')
    n_gm = size(df_gm_info,1);
else
    %select ground motion info file
    msg_inp=sprintf('Select ground motion file [y/N]: '); f_inp = lower(input(msg_inp,'s'));
    switch f_inp
        case 'y'
            [f,p] = uigetfile('*.csv','Select Ground Motion Processing File');
            fn_gm_info = [p,f];
        case 'n'
            fn_gm_info = [fn_eq,'_gm_info'];
            if flag_rot == 1; fn_gm_info = [fn_gm_info,'_rotated']; end
            fn_gm_info = [dir_inpt,fn_gm_info,'.csv'];
    end
    %read time-history info
    df_gm_info       = readtable(fn_gm_info,'VariableNamingRule','preserve');
    %reformat time-history info columns
    if strcmp(class(df_gm_info.eq_id),'double');   df_gm_info.eq_id   = int16(df_gm_info.eq_id);    end
    if strcmp(class(df_gm_info.station),'double'); df_gm_info.station = arrayfun(@(x) num2str(x), df_gm_info.station,'UniformOutput',false); end

    %process only selected time histories
    if any( strcmp('flag_to_process', df_gm_info.Properties.VariableNames) )
        df_gm_info = df_gm_info(df_gm_info.flag_to_process > 0,:);
    end
    %number of time histories to process
    n_gm = size(df_gm_info,1);

    %initialize summary structures 
    df_gm_prc  = table('Size',[n_gm*length(cmp),0]); 
    gm_prc_all  = cell(n_gm,length(cmp));
    gm_seed_all = cell(n_gm,length(cmp));
    gm_corr_all = cell(n_gm,length(cmp));
end

%initialize output files
save([dir_out,fn_gm_prcd],'df_gm_info','df_gm_prc','gm_prc_all','gm_seed_all','gm_corr_all')


%% Processing
%initialize gm index
k = 0;
%iterate over ground motions
while true
    fprintf('Processed ground motion %i of %i\n',k+1,n_gm)
    while true
        msg_inp=sprintf('Select ground motion for processing [next:n, specific: gm index, terminate: e]: '); f_prc = lower(input(msg_inp,'s'));
        switch f_prc
            case 'e';  k = -1;
            case 'n';  k = k+1;
            otherwise; k = str2double(f_prc);
        end
        %confirm selection
        if     strcmp(f_prc,'n') && k <=n_gm %move to the next ground motion
            flag_cont = 'y';
        elseif strcmp(f_prc,'n') && k >n_gm  %exit, no ground motions left for processing
            flag_cont = 'y'; k =-1;
        elseif k >0 && k <=n_gm              %move to sepecific ground motion
            flag_cont = input(sprintf('Process %s-%s [y/N]: ',df_gm_info{k,'network'}{1},df_gm_info{k,'station'}{1}),'s');
        elseif k == -1                       %exit after confirmation
            flag_cont = input('Exit [y/N]: ','s');
        else                                 %invalid input
            flag_cont = 'n';
        end
        %exit loop
        if strcmpi(flag_cont,'y'); break; end
    end
    %terminate processing
    if k == -1; break; end

    %ground motion name
    fn_gm = df_gm_info{k,sprintf('fname_%s',cmp{1})}{1};
    fn_gm = fn_gm(1:strfind(fn_gm,['_',cmp{1}])-1);

    %load time histories
    fn_gm_rec = cell(1,length(cmp));
    data      = cell(1,length(cmp));
    gm_raw    = cell(1,length(cmp));
    for j = 1:length(cmp)
        %create file name
        fn_gm_rec{j} = df_gm_info{k,sprintf('fname_%s',cmp{j})}{1};
        %read time history
        data{j} = readtable([dir_inpt,fn_gm_rec{j}],'FileType','text');
        %compute time history
        gm_raw{j} = calc_his(data{j}{:,:},per,'acc',acc_grav);
    end
    %ground motion info
    n_sta = sprintf('%s-%s',df_gm_info{k,'network'}{1}, df_gm_info{k,'station'}{1});
    r_hyp = df_gm_info{k,'hyp_dist'};
    r_rup = nan;
    if ismember('r_rup',df_gm_info.Properties.VariableNames); r_rup = df_gm_info{k,'r_rup'}; end

    %baseline correction
    flag_exit = 'reset';
    while ~strcmpi(flag_exit,'y')
        %choose input method
        while true
            flag_inp = input('Input Method (Keyboard=1, Cursor=2): ');
            if ismember(flag_inp,[1,2]); break; end
        end

        %pre-event correction
        t_pre = 0;
        acc_off = zeros(1,length(cmp));

        %baseline correction loop
        while ~strcmpi(flag_exit,'y')
            %initialize
            t1 = nan; t2 = nan; t3 = nan; ts = nan; 
            tsearch = nan; fc = nan; 
            maxE = nan;
            K = nan(1,length(cmp)); Km = nan(1,length(cmp)); 
            t1f = nan(1,length(cmp)); t2f = nan(1,length(cmp));
            %inpupt/processing
            if strcmpi(flag_exit,'reset')
                gm_prc  = cell(1,length(cmp));
                gm_corr = cell(1,length(cmp));
                gm_seed = cell(1,length(cmp));
                text_inp = ' ';
                i_attempt = -2;
            elseif ~strcmpi(flag_exit,'y')
                %intialize ground motion (truncation and detrending)
                if i_attempt == -1
                    %truncate motion information
                    flag_trc = input('truncate ground motion [y/N]: ','s'); 
                    t_trnc_s=0; t_trnc_e= inf;
                    if strcmpi(flag_trc,'y')
                        msg_inp='select start point (sec): '; t_trnc_s = input(msg_inp);
                        msg_inp='select end point   (sec): '; t_trnc_e = input(msg_inp);
                    end
                    %truncate time histories
                    i_t = and(data{1}{:,1}>=t_trnc_s, data{1}{:,1}<=t_trnc_e);
                    time_array = data{1}{i_t,1};
                    acc_array = nan(length(time_array),3);
                    for j = 1:length(cmp)
                        acc_array(:,j) = data{j}{i_t,2};
                    end
                    %initialize processed gm 
                    acc_prc  = zeros(length(time_array),length(cmp)); 
                    acc_corr = zeros(length(time_array),length(cmp)); 
                elseif i_attempt == 0
                    %compute acc offset
                    flag_acc = input('remove acc offset [y/N]: ','s'); 
                    acc_off = zeros(1,length(cmp));
                    if strcmpi(flag_acc,'y')
                        msg_inp='select pre-event portion (sec): '; if flag_inp==1; t_pre = input(msg_inp); elseif flag_inp==2; fprintf('%s\n',msg_inp); [t_pre, ~] = ginput(1); end
                        for j = 1:length(cmp)
                            acc_off(j) = mean( acc_array(time_array<=t_pre,j) );
                        end
                    end
                    %remove acc offset
                    for j = 1:length(cmp)
                        acc_array(:,j) = acc_array(:,j) - acc_off(j);
                    end                
                end

                %baseline correction method
                if i_attempt < 1
                    flag_method = 0;
                    %seed ground motion
                    for j = 1:length(cmp)
                        gm_seed{j} = calc_his([time_array,acc_array(:,j)],per,'acc',acc_grav);
                    end
                else
                    flag_method = input(sprintf(['Baseline correction options:\n',...
                                                 '\t1: Boore Piecewise      \n\t2: Boore Quad \n', ...
                                                 '\t3: High-pass Filter     \n\t4: Static     \n', ...
                                                 '\t5: Static Search:       \n\t6: Step Fit   \n',...
                                                 '\t7: Step Fit Constrained \n\t8: Trilinear  \n', ...
                                                 '  Select method: ']));
                end
                
                switch flag_method
                    case 0
                        method = 'Initialization';
                        for j = 1:length(cmp) 
                            acc_prc(:,j) = acc_array(:,j);
                        end
                    case 1 % [d v a]=baselineK(t,ai,'BoorePieceWise',t1,t2)
                        %baseline correction method 
                        method = 'BoorePieceWise';
                        %collect input
                        msg_inp='baseline correction start time (sec): ';        if flag_inp==1; t1 = input(msg_inp); elseif flag_inp==2; fprintf('%s\n',msg_inp); [t1, ~] = ginput(1); end
                        msg_inp='baseline correction intermediate time (sec): '; if flag_inp==1; t2 = input(msg_inp); elseif flag_inp==2; fprintf('%s\n',msg_inp); [t2, ~] = ginput(1); end
                        %input text
                        text_inp = sprintf('Method: %s, t1: %.1f, t2: %.1f',method,t1,t2);
                        %apply baseline correction
                        for j = 1:length(cmp) 
                            [~, ~, acc_prc(:,j), ~, ~, acc_corr(:,j)] = baselineK(time_array-t_trnc_s, acc_array(:,j), method ,t1-t_trnc_s, t2-t_trnc_s);
                        end
                    case 2 % [d v a]=baselineK(t,ai,'BooreQuad',t1)
                        %baseline correction method 
                        method = 'BooreQuad';
                        %collect input
                        msg_inp='baseline correction start time (sec): '; if flag_inp==1; t1 = input(msg_inp); elseif flag_inp==2; fprintf('%s\n',msg_inp); [t1, ~] = ginput(1); end
                        %input text
                        text_inp = sprintf('Method: %s, t1: %.1f',method,t1);
                        %apply baseline correction
                        for j = 1:length(cmp) 
                            [~, ~, acc_prc(:,j), ~, ~, acc_corr(:,j)] = baselineK(time_array-t_trnc_s, acc_array(:,j), method ,t1-t_trnc_s);
                        end
                    case 3 % [d v a]=baselineK(t,ai,'Filter',Fc)
                        %baseline correction method 
                        method = 'Filter';
                        %collect input
                        msg_inp = 'baseline high-pass corner frequency (hz): '; fc = input(msg_inp);
                        %input text
                        text_inp = sprintf('Method: %s, fc: %.2f',fc);
                        %apply baseline correction
                        for j = 1:length(cmp) 
                            [~, ~, acc_prc(:,j), ~, ~, acc_corr(:,j)] = baselineK(time_array-t_trnc_s, acc_array(:,j), method ,fc);
                        end
                    case 4 % [d v a t2fit Km]=baselineK(t,ai,'Static',t1,t2,K)
                       %baseline correction method 
                        method = 'Static';
                        %collect input
                        msg_inp='baseline correction start time (sec): ';        if flag_inp==1; t1 = input(msg_inp); elseif flag_inp==2; fprintf('%s\n',msg_inp); [t1, ~] = ginput(1); end
                        msg_inp='baseline correction intermediate time (sec): '; if flag_inp==1; t2 = input(msg_inp); elseif flag_inp==2; fprintf('%s\n',msg_inp); [t2, ~] = ginput(1); end
                        for j = 1:length(cmp) 
                            msg_inp=sprintf('static field constraint, cmp: %s (m): ',cmp{j}); K(j) = input(msg_inp);
                        end
                        %input text
                        text_inp = sprintf('Method: %s, t1: %.1f, t2: %.1f, K: %s',method,t1,t2,num2str(K));
                        %apply baseline correction
                        for j = 1:length(cmp) 
                            [~, ~, acc_prc(:,j), ~, ~, acc_corr(:,j), t2f(j), Km(j)] = baselineK(time_array-t_trnc_s, acc_array(:,j), method ,t1-t_trnc_s, t2-t_trnc_s, K);
                        end
                    case 5 % [d v a t2fit Km]=baselineK(t,ai,'StaticSearch',t1,t2,tsearch,K)
                       %baseline correction method 
                        method = 'StaticSearch';
                        %collect input
                        msg_inp='baseline correction start time (sec): ';        if flag_inp==1; t1 = input(msg_inp); elseif flag_inp==2; fprintf('%s\n',msg_inp); [t1, ~] = ginput(1); end
                        msg_inp='baseline correction intermediate time (sec): '; if flag_inp==1; t2 = input(msg_inp); elseif flag_inp==2; fprintf('%s\n',msg_inp); [t2, ~] = ginput(1); end
                        msg_inp='time for averaging in search (sec): ';          tsearch = input(msg_inp);
                        for j = 1:length(cmp) 
                            msg_inp=sprintf('static field constraint, cmp: %s (m): ',cmp{j}); K(j) = input(msg_inp);
                        end
                        %input text
                        text_inp = sprintf('Method: %s, t1: %.1f, t2: %.1f, tsearch: %.1f, K: [%s]',method,t1,t2,tsearch,num2str(K));
                        %apply baseline correction
                        for j = 1:length(cmp) 
                            [~, ~, acc_prc(:,j), ~, ~, acc_corr(:,j), t2f(j), Km(j)] = baselineK(time_array-t_trnc_s, acc_array(:,j), method ,t1-t_trnc_s, t2-t_trnc_s, tsearch, K);
                        end
                    case 6 % [d v a t1fit t2fit]=baselineK(t,ai,'StepFit',tstart,maxE)
                       %baseline correction method 
                        method = 'StepFit';
                        %collect input
                        msg_inp='start of motion, usually a p-wave pick (sec): ';  if flag_inp==1; ts = input(msg_inp); elseif flag_inp==2; fprintf('%s\n',msg_inp); [ts, ~] = ginput(1); end
                        msg_inp='proportion of energy accumulated [0,1]: ';        maxE = input(msg_inp);
                        %input text
                        text_inp = sprintf('Method: %s, ts: %.1f, maxE: %.2f',method,ts,maxE);
                        %apply baseline correction
                        for j = 1:length(cmp) 
                            [~, ~, acc_prc(:,j), ~, ~, acc_corr(:,j), t1f(j), t2f(j)] = baselineK(time_array-t_trnc_s, acc_array(:,j), method , ts-t_trnc_s, maxE);
                        end
                    case 7 % [d v a t1fit t2fit]=baselineK(t,ai,'StepFit',tstart,maxE,K)
                       %baseline correction method 
                        method = 'StepFit';
                        %collect input
                        msg_inp='start of motion, usually a p-wave pick (sec): '; if flag_inp==1; ts = input(msg_inp); elseif flag_inp==2; [ts, ~] = ginput(1); end
                        msg_inp='proportion of energy accumulated [0,1]: ';       maxE = input(msg_inp);
                        for j = 1:length(cmp) 
                            msg_inp=sprintf('static field constraint, cmp: %s (m): ',cmp{j}); K(j) = input(msg_inp);
                        end
                        %input text
                        text_inp = sprintf('Method: %s, ts: %.1f, maxE: %.2f, K: %.2f',method,ts,maxE,K);
                        %apply baseline correction
                        for j = 1:length(cmp)
                            [~, ~, acc_prc(:,j), ~, ~, acc_corr(:,j), t1f(j), t2f(j)] = baselineK(time_array-t_trnc_s, acc_array(:,j), method , ts-t_trnc_s, maxE, K);
                        end
                    case 8 % trilinear_baseline_fit(t,a_inp,t1,t2,t3)
                        %baseline correction method 
                        method = 'Trilinear';
                        %collect input
                        msg_inp='baseline correction start time (sec): ';        if flag_inp==1; t1 = input(msg_inp); elseif flag_inp==2; fprintf('%s\n',msg_inp); [t1, ~] = ginput(1); end
                        msg_inp='baseline correction intermediate time (sec): '; if flag_inp==1; t2 = input(msg_inp); elseif flag_inp==2; fprintf('%s\n',msg_inp); [t2, ~] = ginput(1); end
                        msg_inp='baseline correction finel time (sec): ';        if flag_inp==1; t3 = input(msg_inp); elseif flag_inp==2; fprintf('%s\n',msg_inp); [t3, ~] = ginput(1); end
                        %input text
                        text_inp = sprintf('Method: %s, t1: %.1f, t2: %.1f, t3: %.1f',method,t1,t2,t3);
                        %apply baseline correction
                        for j = 1:length(cmp) 
                            [~, ~, acc_prc(:,j), ~, ~, acc_corr(:,j)] = trilinear_baseline_fit(time_array-t_trnc_s, acc_array(:,j), t1-t_trnc_s,t2-t_trnc_s,t3-t_trnc_s);
                        end
                    otherwise
                        continue
                end

                if i_attempt > 0
                    %fiter ground motions
                    flag_taper = input('filter ground motion [y/N]: ','s');
                    if strcmpi(flag_taper,'y')
                        %collect input
                        msg_inp = 'baseline high-pass corner frequency (hz): '; fc = input(msg_inp);
                        %filter ground motions
                        for j = 1:length(cmp) 
                            [~, ~, acc_prc(:,j), ~, ~, a_c] = baselineK(time_array-t_trnc_s, acc_prc(:,j), 'Filter' ,fc);
                            acc_corr(:,j) = acc_corr(:,j) + a_c;
                        end
                    end
    
                    %apply tapering
                    flag_taper = input('taper ground motion [y/N]: ','s');
                    if strcmpi(flag_taper,'y')
                        %collect input
                        msg_inp = 'taper ratio in the beginning: '; tpr_str = input(msg_inp);
                        msg_inp = 'taper ratio in the end: ';       tpr_end = input(msg_inp);
                        %tapering
                        taper =  taper_halfsin(time_array-t_trnc_s,tpr_str,tpr_end);
                        for j = 1:length(cmp); acc_prc(:,j) = taper .* acc_prc(:,j); end
                    end
                end

                %compute attributes of baseline corrected time history
                gm_prc  = cell(length(cmp),1);
                gm_corr = cell(length(cmp),1);
                for j = 1:length(cmp)
                    %general ground motion information
                    gm_prc{j} = struct('event', df_gm_info{k,'event'}{1}, 'eq_id',df_gm_info{k,'eq_id'}, 'mag',df_gm_info{k,'mag'}, ...
                                       'network', df_gm_info{k,'network'}{1}, 'station', df_gm_info{k,'station'}{1}, ...
                                       'file_name', fn_gm_rec{j}, 'component',cmp{j});
                    %time history information 
                    gm_prc_partial = calc_his([time_array,acc_prc(:,j)],per,'acc',acc_grav);
                    %copy gm_prc_partial fields to gm_prc
                    f = fieldnames(gm_prc_partial);
                    for jj = 1:length(f); gm_prc{j}.(f{jj}) = gm_prc_partial.(f{jj}); end
                    %add truncation infomation
                    gm_prc{j}.t_trnc_start = t_trnc_s;
                    gm_prc{j}.t_trnc_end   = t_trnc_e;
                    %add offset correction info
                    gm_prc{j}.tp = t_pre-t_trnc_s;
                    gm_prc{j}.acc_offset = acc_off(j);
                    %add baseline correction information
                    gm_prc{j}.method = method;
                    gm_prc{j}.t1 = t1-t_trnc_s; gm_prc{j}.t2 = t2-t_trnc_s; gm_prc{j}.t3 = t3-t_trnc_s; 
                    gm_prc{j}.ts = ts-t_trnc_s; gm_prc{j}.tsearch = tsearch-t_trnc_s; 
                    gm_prc{j}.fc = fc; gm_prc{j}.maxE = maxE;
                    gm_prc{j}.K = K(j);
                    gm_prc{j}.t1f = t1f(j); gm_prc{j}.t2f = t2f(j); gm_prc{j}.Km = Km(j);
                    %corretion time history
                    % gm_corr{j} = calc_his([time_array,acc_corr(:,j)-acc_off(j)],per,'acc',acc_grav);
                    gm_corr{j} = calc_his([time_array,acc_corr(:,j)],per,'acc',acc_grav);
                end
            end

            %plot title
            if isnan(r_rup); text_title = {sprintf('Ground-motion: %s [%i of %i] (R_{hyp}=%.1fkm)',                 strrep(fn_gm,'_','\_'),k,n_gm,r_hyp);       text_inp};
            else;            text_title = {sprintf('Ground-motion: %s [%i of %i] (R_{hyp}=%.1fkm, R_{rup}=%.1fkm)', strrep(fn_gm,'_','\_'),k,n_gm,r_hyp,r_rup); text_inp};
            end
            %produce group plots
            if strcmpi(flag_exit,'reset')
                fig_main = plot_gm2(text_title,gm_raw,[],cmp,acc_unt,{'Seed Motion'},     {'#D95319'},10);
                flag_exit = 'n';
            elseif ~strcmpi(flag_exit,'y')
                fig_main = plot_gm2(text_title,gm_prc,[],cmp,acc_unt,{'Processed Motion'},{'#0072BD'},10);
            end

            %exit flag
            if i_attempt > 0
                flag_exit = input('exit current ground motion [y/N/reset/compound]: ','s'); 
            end
            i_attempt = i_attempt + 1;

            %rest and repeate options
            if strcmpi(flag_exit,'reset') 
                fprintf(' iterating component\n'); 
                break;
            elseif strcmpi(flag_exit,'n') && i_attempt > 0
                fig_iter = plot_gm2(text_title,gm_seed,gm_prc,cmp,acc_unt,{'Seed Motion','Processed Motion'},{'#D95319','#0072BD'},11);
            elseif strcmpi(flag_exit,'compound') && i_attempt > 0
                disp("compounding")
                fig_iter = plot_gm2(text_title,gm_seed,gm_prc,cmp,acc_unt,{'Seed Motion','Updated Seed Motion'},{'#D95319','#0072BD'},11);
                gm_seed = gm_prc; %update seed for compound correction
                for j = 1:length(cmp)
                    acc_array(:,j) = gm_seed{j}.acc;
                end
                flag_exit = 'n';
            end

        end
    end

    %store gm structure
    for j = 1:length(cmp)
        gm_prc_all{k,j}  = gm_prc{j};
        gm_seed_all{k,j} = gm_seed{j};
        gm_corr_all{k,j} = gm_corr{j};
    end

    %summarize ground motion
    for j = 1:length(cmp)
        kk = 3*(k-1) + j;
        %general info
        cn_gm_info = {'eq_id','event','mag','hyp_lat','hyp_long','hyp_depth','network','station','sta_lat','sta_long','sta_elev','hyp_dist'};
        df_gm_prc(kk,cn_gm_info) = df_gm_info(k, cn_gm_info);
        %gm information
        df_gm_prc{kk,'component'}         = cmp(j);
        df_gm_prc{kk,{'dt','npt'}}        = [gm_prc{j}.dt,  gm_prc{j}.npt];
        df_gm_prc{kk,{'PGA','PGV','PGD'}} = [gm_prc{j}.pga, gm_prc{j}.pgv, gm_prc{j}.pgd];
        df_gm_prc{kk,'disp_final1'}       = gm_prc{j}.dis(end);
        df_gm_prc{kk,'disp_final2'}       = mean( gm_prc{j}.dis(gm_prc{j}.time > gm_prc{j}.time(end)-1) );
        %processing information
        df_gm_prc{kk,'method'}                       = {method};
        df_gm_prc{kk,{'tp'}}                         = t_pre-t_trnc_s;
        df_gm_prc{kk,{'t1','t2','t3','tstart','tsearch'}} = [t1,t2,t3,ts,tsearch]-t_trnc_s;
        df_gm_prc{kk,{'HP-FC','K','maxE'}}           = [fc,K(j),maxE];
        df_gm_prc{kk,{'t1f','t2f','Km'}}             = [t1f(j),t2f(j),Km(j)];
        %extract raw ground motion time
        gm_time_raw   = regexp(fn_gm,'[0-9]+','match','once');
        %sufix gm name
        fn_gm_sufix = fn_gm( strfind(fn_gm,gm_time_raw)+length(gm_time_raw) : end );
        %year, month, day
        gm_time_year  = str2double(gm_time_raw(1:4));
        gm_time_month = str2double(gm_time_raw(5:6));
        gm_time_day   = str2double(gm_time_raw(7:8));
        %hour, minutes, sec
        gm_time_hour = str2double(gm_time_raw(9:10));
        gm_time_min  = str2double(gm_time_raw(11:12));
        gm_time_sec  = str2double(gm_time_raw(13:14));
        %offset gm time based on start truncation
        gm_time_sec  = gm_time_sec + t_trnc_s;
        %fix sec [0,60] carry over to min
        gm_time_min  = gm_time_min  + floor(gm_time_sec/60);
        gm_time_sec  = mod(gm_time_sec, 60);
        %fix min [0,60] carry over to hour
        gm_time_hour = gm_time_hour + floor(gm_time_min/60);
        gm_time_min  = mod(gm_time_min, 60);
        %fix hour [0,24] carry over to day
        gm_time_day  = gm_time_day  + floor(gm_time_hour/60);
        gm_time_hour = mod(gm_time_hour, 24);
        %general filename
        fn_gm_rec = sprintf('%.4i%.2i%.2i%.2i%.2i%.2i%s_%s',gm_time_year,gm_time_month,gm_time_day,gm_time_hour,gm_time_min,gm_time_sec, fn_gm_sufix,cmp{j});
        %gm time array
        time = arrayfun(@(x) num2str(x,'%012.8f'), gm_prc{j}.time-t_trnc_s, 'UniformOutput', 0);
        %acceleration
        df_gm_prc{kk,{'fname_acc'}} = {[fn_gm_rec,'.acc']};
        acc  = arrayfun(@(x) num2str(x,'% 13f'), gm_prc{j}.acc, 'UniformOutput', 0);
        acc = cell2table([time,acc],'VariableNames',{'time','acc'});
        writetable(acc,[dir_out,df_gm_prc{kk,{'fname_acc'}}{1}],'Delimiter',' ','QuoteStrings','none','WriteVariableNames',false,'FileType','text')
        %velocity
        df_gm_prc{kk,{'fname_vel'}} = {[fn_gm_rec,'.vel']};
        vel  = arrayfun(@(x) num2str(x,'% 13f'), gm_prc{j}.vel, 'UniformOutput', 0);
        vel = cell2table([time,vel],'VariableNames',{'time','vel'});
        writetable(vel,[dir_out,df_gm_prc{kk,{'fname_vel'}}{1}],'Delimiter',' ','QuoteStrings','none','WriteVariableNames',false,'FileType','text')
        %displacement
        df_gm_prc{kk,{'fname_disp'}} = {[fn_gm_rec,'.disp']};
        disp  = arrayfun(@(x) num2str(x,'% 13f'), gm_prc{j}.dis, 'UniformOutput', 0);
        disp = cell2table([time,disp],'VariableNames',{'time','dis'});
        writetable(disp,[dir_out,df_gm_prc{kk,{'fname_disp'}}{1}],'Delimiter',' ','QuoteStrings','none','WriteVariableNames',false,'FileType','text')

        %produce individual plots
        %plot title
        if isnan(r_rup); text_title = {sprintf('Seed Motion: %s (R_{hyp}=%.1fkm)',                 strrep(fn_gm_rec,'_','\_'), r_hyp);       text_inp};
        else;            text_title = {sprintf('Seed Motion: %s (R_{hyp}=%.1fkm, R_{rup}=%.1fkm)', strrep(fn_gm_rec,'_','\_'), r_hyp,r_rup); text_inp};
        end
        %plot original and processed time histories
        [fig1, fig2] = plot_gm(text_title,gm_seed{j},gm_prc{j},gm_corr{j},'m/sec^2');
        [fig3, fig4] = plot_gm(text_title,gm_prc{j},[],[],                'm/sec^2',{'Processed Motion'},{'#0072BD'},3);
        %save figures
        fn_fig  = [dir_fig,fn_gm_rec,'.pdf'];
        fn_fig1 = [dir_fig,fn_gm_rec,'1.pdf']; saveas(fig1,fn_fig1)
        fn_fig2 = [dir_fig,fn_gm_rec,'2.pdf']; saveas(fig2,fn_fig2)
        fn_fig3 = [dir_fig,fn_gm_rec,'3.pdf']; saveas(fig3,fn_fig3)
        fn_fig4 = [dir_fig,fn_gm_rec,'4.pdf']; saveas(fig4,fn_fig4)
        %merging through pdftk
        %system(sprintf('pdftk %s %s %s %s cat output %s',fn_fig1,fn_fig2,fn_fig3,fn_fig4,fn_fig));
        %merging through matlab
        mergepdf({fn_fig1,fn_fig2,fn_fig3,fn_fig4},fn_fig)
        system(sprintf('rm %s %s %s %s',fn_fig1,fn_fig2,fn_fig3,fn_fig4));
    end
    %save main figure
    fn_gm = sprintf('%.4i%.2i%.2i%.2i%.2i%.2i%s_all',gm_time_year,gm_time_month,gm_time_day,gm_time_hour,gm_time_min,gm_time_sec, fn_gm_sufix);   
    fn_fig_main = [dir_fig,fn_gm,'.pdf']; saveas(fig_main,fn_fig_main)

    %save intermediate results
    save([dir_out,fn_gm_prcd],'df_gm_info','df_gm_prc','gm_prc_all','gm_seed_all','gm_corr_all')
end

%% Output
%keep only processed ground motions
switch class(df_gm_prc.eq_id)
    case 'double'
        df_gm_prc = df_gm_prc(df_gm_prc.eq_id<0,:);
    case 'cell'
        df_gm_prc = df_gm_prc(~cellfun(@(x) isempty(x), df_gm_prc.eq_id),:);
end
%save gm processing info
fn_gm_info = [fn_eq,'_gm_info_tilt_corrected.csv'];
writetable(df_gm_prc,[dir_out,fn_gm_info])
%save ground motions
save([dir_out,fn_gm_prcd],'df_gm_info','df_gm_prc','gm_prc_all','gm_seed_all','gm_corr_all')
