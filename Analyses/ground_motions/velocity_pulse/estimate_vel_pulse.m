%% Estimate Velocity Pules 
%
%   This script estimates the parameters of velocity pulses in the baseline
%   corrected time histories

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

%acceleration uits
acc_grav = 9.80665;
acc_unt  = 'm/sec^2';

%components
if     flag_rot == 0; cmp = {'Z','N','E'};
elseif flag_rot == 1; cmp = {'Z','FN','FP'};
end

%psa periods
per = logspace(-2,1,300);

%importance factor disp to vel
fact_imp_d2v = 2;

%theta:    [f_p,  nu, gamma]
theta_0  = [0.5,   0.,  1.0];
theta_lb = [0.05, -pi, 0.25];
theta_ub = [5.0,   pi, 5.0];

% directories
% ---   ---   ---   ---   ---
%input files and directories
dir_inpt = '../../../Data/ground_motions/';
if     flag_rot == 0; dir_inpt = [dir_inpt,'corrected_gm/'];
elseif flag_rot == 1; dir_inpt = [dir_inpt,'corrected_gm_rot/'];
end
if     flag_io == 1
    fn_eq    = '2022_Guanshan';
    dir_inpt = [dir_inpt,'M6.5_0917/'];
elseif flag_io == 2
    fn_eq    = '2022_Chihshang';
    dir_inpt = [dir_inpt,'M6.9_0918/'];
end

%output directories
dir_out = '../../../Data/ground_motions/';
if     flag_rot == 0; dir_out = [dir_out,'vel_pulse/'];
elseif flag_rot == 1; dir_out = [dir_out,'vel_pulse_rot/'];
end
if     flag_io == 1; dir_out = [dir_out,'M6.5_0917/'];
elseif flag_io == 2; dir_out = [dir_out,'M6.9_0918/'];
end
dir_fig = [dir_out,'figures/'];

%ground-motion component information
fn_gm_cmp_info = [fn_eq,'_gm_info_tilt_corrected.csv'];

%optimization parametets
% ---   ---   ---   ---   ---
n_opt_iter = 3;
fmcopt = optimoptions('fmincon','Display','iter','Algorithm','sqp', ...
                      'MaxFunctionEvaluations',5000, 'MaxIterations', 2500, ...
                      'ConstraintTolerance', 1.0000e-09,'StepTolerance', 1.0000e-12);

%% Read Input
%read time-history info
df_th_info       = readtable([dir_inpt,fn_gm_cmp_info],'VariableNamingRule','preserve');
%reformat time history info table
if strcmp(class(df_th_info.eq_id),'double');   df_th_info.eq_id   = int16(df_th_info.eq_id);    end
if strcmp(class(df_th_info.station),'double'); df_th_info.station = arrayfun(@(x) num2str(x), df_th_info.station,'UniformOutput',false); end
n_th = size(df_th_info,1);

%ground motion info columns
cn_gm_info = {'eq_id','event','mag','hyp_lat','hyp_long','hyp_depth','network','station','sta_lat','sta_long','sta_elev','hyp_dist'};

%ground motion info
[~,i_gm] = unique(df_th_info(:,cn_gm_info(1:end-4)),'rows');
df_gm_info = df_th_info(i_gm,cn_gm_info);
n_gm = size(df_gm_info,1);
assert(n_gm*length(cmp)==n_th,'Error. Inconsistent ground motion and time history number.')

%short time-history info for consistency
[~,i_th] = sortrows(df_th_info(:,cn_gm_info(1:end-4)) );
df_th_info = df_th_info(i_th,:);

%% Processing
%create directories
if not(isfolder(dir_out)); mkdir(dir_out); end
if not(isfolder(dir_fig)); mkdir(dir_fig); end

%initialize summary structures 
df_vel_pulse = df_th_info(:,cn_gm_info); 
gm_inst_all = cell(n_gm,length(cmp));
vel_pulse_all  = cell(n_gm,length(cmp));

%initialize columns
% df_vel_pulse{:,{'dt','npt','PGA','PGV','PGD','disp_res1','disp_res2'}} = nan;
% df_vel_pulse{:,{'t0','A','fp','nu','gamma'}}                           = nan;
% df_vel_pulse{:,{'vp_disp_res','vp_disp_res_ratio','vp_disp_res_sign'}} = nan;
% df_vel_pulse{:,{'fname_acc','fname_vel','fname_disp'}} = repmat({''},n_th,3);
df_vel_pulse{:,{'cmp'}}                                 = {''};
df_vel_pulse{:,{'t0','A','fp','nu','gamma'}}            = nan;
df_vel_pulse{:,{'IvNorm'}}                              = nan;
df_vel_pulse{:,{'vp_vpeak','vp_vp2p'}}                  = nan;
df_vel_pulse{:,{'vp_dpeak','vp_dres','vp_dres_ratio'}}  = nan;
df_vel_pulse{:,{'fname_acc','fname_vel','fname_disp'}}  = repmat({''},n_th,3);


%iterate over ground motions
for k = 1:n_gm
    %ground motion information
    eq_id  = df_gm_info{k,'eq_id'};
    eq_n   = df_gm_info{k,'event'};
    net_n  = df_gm_info{k,'network'}{1};
    sta_n  = df_gm_info{k,'station'}{1};
    nsta_n = sprintf('%s-%s',net_n, sta_n);
    r_hyp  = df_gm_info{k,'hyp_dist'};
    r_rup = nan;
    if ismember('r_rup',df_gm_info.Properties.VariableNames); r_rup = df_gm_info{k,'r_rup'}; end
    %report progress
    fprintf('Processing GM: %s (%i of %i)  ... \n',nsta_n,k,n_gm)

    %ground motion indices
    i_gm = ismember(df_th_info(:,{'eq_id','network','station'}), df_gm_info(k,{'eq_id','network','station'}), 'rows');

    %initalize time history placeholders
    gm_inst   = cell(length(cmp),1);
    vel_pulse = cell(length(cmp),1);

    %iterate over components
    for j = 1:length(cmp)
        %component 
        fprintf(' Component: %s\n',cmp{j})

        %time history index
        i_th = and(i_gm, strcmp(df_th_info.component,cmp{j}));
        th_info = df_th_info(i_th,:);

        %time history name
        fn_gm  = th_info.fname_acc{1}(1:max(strfind(th_info.fname_acc{1},['_',cmp{j}]))-1);
        %read time history
        acc_th = table2array( readtable([dir_inpt, th_info.fname_acc{1}],'FileType','text') );
        %compute intesity parameters
        gm_inst{j} = calc_his(acc_th,per,'acc',acc_grav);

        %search window
        switch th_info.method{1}
            case 'BoorePieceWise'
                t_str = th_info.t1;
                t_end = th_info.t2;
            case 'BooreQuad'
                t_str = th_info.t1;
                t_end = acc_th(end,1);
            case 'Filter'
                error('Error: Unimplemented method')
            case 'Static'
                error('Error: Unimplemented method')
            case 'StaticSearch'
                error('Error: Unimplemented method')
            case 'StepFit'
                error('Error: Unimplemented method')
            case 'Trilinear'
                t_str = th_info.t1;
                t_end = th_info.t3;
        end

        %objective function
        fun_obj = @(theta) DetermineAmpTime(gm_inst{j}.time,gm_inst{j}.vel,gm_inst{j}.dis,theta(1),theta(2),theta(3),t_str,t_end,fact_imp_d2v);

        %constrained optimization
        for jj = 1:n_opt_iter
            fprintf(' Optimization Iteration: %i of %i\n',jj,n_opt_iter)
            theta_0(jj+1,:) = fmincon(fun_obj,theta_0(jj,:),[],[],[],[],theta_lb,theta_ub,[],fmcopt);
        end
        theta = theta_0(end,:);

        %velocity pulse parameters
        [mse,A,t0] = fun_obj(theta);
        fp = theta(1); nu = theta(2); gamma = theta(3);

        %phase wrapping for neg amplitude
        if A < 0
            A  = abs(A);
            nu = mod(pi + nu, 2*pi);
        end

        %estimated velocity pulse
        vel_p        = VelPulseMP03(gm_inst{j}.time,A,fp,nu,gamma,t0);
        vel_pulse{j} = calc_his([gm_inst{j}.time,vel_p],per,'vel',acc_grav);

        %identify pulse time
        [~,idx_t] = min(abs(gm_inst{j}.time - t0));

        %time history component
        df_vel_pulse{i_th,{'cmp'}}                      = cmp(j);
        %velocity pulse information
        df_vel_pulse{i_th,{'t0','A','fp','nu','gamma'}} = [t0,A,fp,nu,gamma];
        %normalized arias velocity
        df_vel_pulse{i_th,{'IvNorm'}}                   = gm_inst{j}.IvNorm_his(idx_t);
        %velocity pulse
        df_vel_pulse{i_th,'vp_vpeak'}                   = max(abs(vel_p));
        df_vel_pulse{i_th,'vp_vp2p'}                    = max(vel_p) - min(vel_p);
        %permanent offset
        df_vel_pulse{i_th,'vp_dpeak'}                   = max(abs(vel_pulse{j}.dis));
        df_vel_pulse{i_th,'vp_dres'}                    = abs(vel_pulse{j}.dis(end));
        %permanent displacement ratio
        df_vel_pulse{i_th,'vp_dres_ratio'}              = df_vel_pulse{i_th,'vp_dres'} / df_vel_pulse{i_th,'vp_dpeak'};

        %store time histories
        %acceleration
        df_vel_pulse{i_th,{'fname_acc'}} = {[fn_gm,'_',cmp{j},'.acc']};
        acc = array2table([gm_inst{j}.time,gm_inst{j}.acc,vel_pulse{j}.acc], ...
                          'VariableNames',{'time','gm_acc','vpulse_acc'});
        writetable(acc,[dir_out,df_vel_pulse{i_th,{'fname_acc'}}{1}],'Delimiter',' ','QuoteStrings','none','WriteVariableNames',false,'FileType','text')
        %velocity
        df_vel_pulse{i_th,{'fname_vel'}} = {[fn_gm,'_',cmp{j},'.vel']};
        vel = array2table([gm_inst{j}.time,gm_inst{j}.vel,vel_pulse{j}.vel], ...
                          'VariableNames',{'time','gm_vel','vpulse_vel'});
        writetable(vel,[dir_out,df_vel_pulse{i_th,{'fname_vel'}}{1}],'Delimiter',' ','QuoteStrings','none','WriteVariableNames',false,'FileType','text')
        %displacement
        df_vel_pulse{i_th,{'fname_disp'}} = {[fn_gm,'_',cmp{j},'.disp']};
        dis = array2table([gm_inst{j}.time,gm_inst{j}.dis,vel_pulse{j}.dis], ...
                          'VariableNames',{'time','gm_dis','vpulse_dis'});
        writetable(dis,[dir_out,df_vel_pulse{i_th,{'fname_disp'}}{1}],'Delimiter',' ','QuoteStrings','none','WriteVariableNames',false,'FileType','text')
    end

    %store processed time histories
    gm_inst_all(k,:)   = gm_inst;
    vel_pulse_all(k,:) = vel_pulse;

    %plot title
    if isnan(r_rup); text_title = {sprintf('Ground-motion: %s (R_{hyp}=%.1fkm)',                 strrep(fn_gm,'_','\_'),r_hyp)};
    else;            text_title = {sprintf('Ground-motion: %s (R_{hyp}=%.1fkm, R_{rup}=%.1fkm)', strrep(fn_gm,'_','\_'),r_hyp,r_rup)};
    end
    %plot velocity pulse fit
    fig_main = plot_gm2(text_title,gm_inst,vel_pulse,cmp,acc_unt,{'Processed Motion','Velocity Pulse'},{'#0072BD','#D95319'},10);
    %filename
    fn_fig_main = [dir_fig,fn_gm,'_all','.pdf']; saveas(fig_main,fn_fig_main)
    close all
end

%% Output
%save gm info
fn_prc = [fn_eq,'_gm_info_vel_pulse.csv'];
writetable(df_vel_pulse,[dir_out,fn_prc])

%save ground motions
fn_gm = [fn_eq,'_gm_info_vel_pulse.mat'];
save([dir_out,fn_gm],'df_gm_info','df_vel_pulse','gm_inst_all','vel_pulse_all')


%% Auxiliary Functions
function [mse, A, t0] = DetermineAmpTime(time, vel, dis, f_p, nu, gamma, ts, te, fact_imp)
%DetermineAmpTime

    %valid start times
    idx_t = find( and(time>=ts, time<=te) );
    
    %initialze error term and amplitude
    mse_array = nan(size(idx_t));
    amp_array = nan(size(idx_t));
    
    %iterate over valid time instances
    parfor j = 1:length(idx_t)
        %velocity pulse
        vel_p = VelPulseMP03(time, 1.0,f_p,nu,gamma,time(idx_t(j)));
        %diplacement pulse
        dis_p = VelPulseMP03Disp(time, 1.0,f_p,nu,gamma,time(idx_t(j)));
    
        %determined amplitude
        amp_array(j) = [vel_p; fact_imp*dis_p] \ [vel; fact_imp*dis];
    
        %compute normalized resiuals 
        mse_array(j) = norm([vel-amp_array(j)*vel_p; fact_imp*(dis-amp_array(j)*dis_p)]);
    end
    
    %min mean square error value and position
    [mse,i_min] = min(mse_array);
    mse_array   = [time(idx_t), mse_array];
    
    %determined parameters
    A  = amp_array( i_min );
    t0 = time( idx_t(i_min) );

end 
