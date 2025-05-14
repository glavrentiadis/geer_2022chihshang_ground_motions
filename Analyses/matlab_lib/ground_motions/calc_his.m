function [gm_obj] = calc_his(his,per,type,grav_acc) 
% for now assume that only thing given is acc and time
% and calculate spectra,vel, and disp

%default input
if nargin < 3; type     = 'acce'; end
if nargin < 4; grav_acc = 1.; end

%constants
gravity = 9.80665;

%gravity correction ratio from acc units to cm/sec^2
grav_ratio = 100 * gravity / grav_acc;

if ~isstruct(his)
    %initialize structure
    gm_obj = struct('npt',length(his),'dt',his(2,1)-his(1,1));
    
    %add period
    gm_obj.per = per;
    
    %compute time history
    switch type
        case 'acc'
            acc       = his(:,2);
            %compute vel & disp convert to cm/sec and cm
            [vel,dis] = acce2veldisp(grav_ratio * acc, gm_obj.dt);
        case 'vel'
            vel        = his(:,2);
            %compute acc and dip
            [acc,dis] = vel2accedisp(vel, gm_obj.dt);
            acc = acc/grav_ratio;
        otherwise
            error('Invalid type option')
    end

    %add time history information
    gm_obj.time = his(:,1);
    gm_obj.acc  = acc;
    gm_obj.vel  = vel;
    gm_obj.dis  = dis;
else
    gm_obj = his;
end

%spectral accelerations
gm_obj.ars = rspec(gm_obj.time,gm_obj.acc,gm_obj.per,0.05)';

%add intensity measures later for plots and other ouput(TODO)
gm_obj.pga = max(abs(gm_obj.acc));
gm_obj.pgv = max(abs(gm_obj.vel));
gm_obj.pgd = max(abs(gm_obj.dis));

% arias intensity
gm_obj = husid(gm_obj,grav_acc);

end