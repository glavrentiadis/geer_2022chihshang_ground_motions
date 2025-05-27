%% Rupture Distance Computation
%
%   This script computes the station distance metrics with respect to the 
%   Central Range and Longitudinal Valley

%add user functions
addpath('../matlab_lib/distances/')

%% Define Input
%input/output flag
%  flag_io = 1 for Central Range Fault
%  flag_io = 2 for Longitudinal Valley Fault
%  flag_io = 3 for Central Range Fault (Rupture Extent)
flag_io = 1;

%fault name
if flag_io == 1;     flt_name = 'CRF';
elseif flag_io == 2; flt_name = 'LVF';
elseif flag_io == 3; flt_name = 'CRF_rupt';
end

%number of discretization points
n_pt = 2500;

%station information
fn_sta_info     = '../../Raw_files/ground_motions/instrument_corrected/station_info.csv';

%finite fault model
if flag_io == 1;     fname_fflt = '../../Raw_files/gis/CRF_finite_fault_model.csv';
elseif flag_io == 2; fname_fflt = '../../Raw_files/gis/LVF_finite_fault_model.csv';
elseif flag_io == 3; fname_fflt = '../../Raw_files/gis/CRF_finite_fault_model_mod.csv';
end
%fault surface projection
if flag_io == 1;     fname_tflt = '../../Raw_files/gis/CRF_finite_fault_trace.csv';
elseif flag_io == 2; fname_tflt = '../../Raw_files/gis/LVF_finite_fault_trace.csv';
elseif flag_io == 3; fname_tflt = '../../Raw_files/gis/CRF_finite_fault_trace_mod.csv';
end

%output directory
dir_out = '../../Data/ground_motions/dist_metrics/';

%% Read Input
%fault trace
df_flttrc = readtable(fname_tflt);
%finite fault model
df_fltfin = readtable(fname_fflt);

%station information
df_sta_info = readtable(fn_sta_info);
n_sta = size(df_sta_info,1);

%% Processing
% coordinate transformation
%------------------------------------------
%project lat/lon coordinates to UTM
%define projection system
utm_zone = '51R';
utmstruct = defaultm('utm'); 
utmstruct.zone = utm_zone;  
utmstruct.geoid = wgs84Ellipsoid;
utmstruct = defaultm(utmstruct);
%calc utm coordinates
%trace fault
fltt_xyz = nan(size(df_flttrc,1),3);
[fltt_xyz,fltt_xyz(:,2)] = projfwd(utmstruct,df_flttrc.latitude,df_flttrc.longitude); %finite Cartesian coordinates
fltt_xyz(:,3)            = 1000 * df_flttrc.depth;
%finte fault
fltf_xyz = nan(size(df_fltfin,1),3);
[fltf_xyz,fltf_xyz(:,2)] = projfwd(utmstruct,df_fltfin.latitude,df_fltfin.longitude); %finite Cartesian coordinates
fltf_xyz(:,3)            = 1000 * df_fltfin.depth;
%station coordinates
sta_xyz = nan(size(df_fltfin,1),3);
[sta_xyz,sta_xyz(:,2)] = projfwd(utmstruct,df_sta_info.latitude,df_sta_info.longitude);
sta_xyz(:,3)            = 0;
%convert to km
fltf_xyz = fltf_xyz/1000;
fltt_xyz = fltt_xyz/1000;
sta_xyz  = sta_xyz/1000;
%fault range
range_xyz = [min(fltf_xyz,[],1); max(fltf_xyz,[],1)];

% create fault plane on dense grid
%------------------------------------------
%create x & z grid
[grid_x,grid_z] = meshgrid(linspace(range_xyz(1,1),range_xyz(2,1),n_pt),linspace(range_xyz(1,3),range_xyz(2,3),n_pt));
grid_x = reshape(grid_x,1,[]);
grid_z = reshape(grid_z,1,[]);
%interpolate y dimension
fun_interp = scatteredInterpolant(fltf_xyz(:,1),fltf_xyz(:,3),fltf_xyz(:,2));
fun_interp.ExtrapolationMethod = 'none';
grid_y = fun_interp(grid_x,grid_z);
%remove unavailable points
i_nan = isnan(grid_y);
grid_x = grid_x(~i_nan);
grid_y = grid_y(~i_nan);
grid_z = grid_z(~i_nan);
%dense rupture grid
fflt_grid_xyz = [grid_x',grid_y',grid_z'];
%dense surface projection grid
fflt_grid_surf_xyz      = fflt_grid_xyz;
fflt_grid_surf_xyz(:,3) = 0;

% compute rupture distance
%------------------------------------------
%rupture distance
r_rup = nan(n_sta,1);
for j = 1:n_sta
    r_rup(j) = min( vecnorm(fflt_grid_xyz-sta_xyz(j,:),2,2) );
end
%joyner-boore distance
r_jb = nan(n_sta,1);
for j = 1:n_sta
    r_jb(j) = min( vecnorm(fflt_grid_xyz(:,1:2)-sta_xyz(j,1:2),2,2) );
end
%GC2 distances
r_x = nan(n_sta,1);
r_y = nan(n_sta,1);
for j = 1:n_sta
    [r_x(j), r_y(j)] = clstdist_point2segments(sta_xyz(j,1:2), fltt_xyz(:,1:2));
end
%dip angle
dip = nan(n_sta,1);
for j = 1:n_sta
    [~,i] = min( vecnorm(fltf_xyz-sta_xyz(j,:),2,2)  );
    dip(j)   = df_fltfin{i,'dip'};
end

% %keep only relevant information
% df_sta_info = df_sta_info(:,{'network','station'});
%store rupture distances
df_sta_info{:,['dip']}   = dip;
df_sta_info{:,['r_rup']} = r_rup;
df_sta_info{:,['r_jb']}  = r_jb;
df_sta_info{:,['r_x']}   = r_x;
df_sta_info{:,['r_y']}   = r_y;


%% Output
%create directories
if not(isfolder(dir_out)); mkdir(dir_out); end

%save rupture distance metrics
fname_sta_info = [flt_name,'_sta_dist_metrics'];
writetable(df_sta_info, [dir_out,fname_sta_info,'.csv']);
