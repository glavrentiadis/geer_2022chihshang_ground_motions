function [taper] = taper_halfsin(time,tpr_str,tpr_end)
%taper_function returns an half sin tapering function
%   tpr_str: percent tapering in the beginning
%   tpr_end: percent tapering in the end

assert(all([tpr_str,tpr_end] <= 1),'Error. Start and end tapers must be below 1')

%start/end times
t_s = min(time);
t_e = max(time);

%tapering frequencies
f_s = 0.5/(tpr_str*(t_e-t_s));
f_e = 0.5/(tpr_end*(t_e-t_s));

%tapering function
taper = ones(length(time),1);
%start tapering
i_tpr = time<=t_s+0.5/f_s;
taper(i_tpr) = 0.5 + 0.5*sin(2*pi*f_s * (t_s + time(i_tpr)) - pi/2);
%end tapering
i_tpr = time>=t_e-0.5/f_e;
taper(i_tpr) = 0.5 + 0.5*sin(2*pi*f_e * (t_e - time(i_tpr)) - pi/2);

end