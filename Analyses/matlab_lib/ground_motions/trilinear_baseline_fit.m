function [d, v, a, d_corr, v_corr, a_corr]= trilinear_baseline_fit(t,a_inp,t1,t2,t3)
%trilinear_baseline_fit performs tirlinear base linecorrection
%important for very long time histories
%
% Input arguments:
%   t:      time array
%   a_inp:  input acceleration time history
%   t1:     intial time
%   t2:     middle time
%   t3:     end time

%convert time and acceleration input to column vectors
t     = reshape(t,    max(size(t)),   min(size(t)));
a_inp = reshape(a_inp, max(size(a_inp)),min(size(a_inp)));
%compute to velocity and displacement of input ground motion
v_inp = cumtrapz(t,a_inp);
d_inp = cumtrapz(t,v_inp);
%end time
tf=max(t);

% Trilinear fit
%---   ---   ---   ---   ---   ---   ---   ---   ---
%determine vel slope at final part of record (t3 - tf)
% . . . . . . . . . . . . . . . . . . . . . . . . . . 
i_fnl = find(t>=t3); %times after t3
%least squares fit
G = [ones(size(v_inp(i_fnl))) t(i_fnl)];
m=lsqlin(G,v_inp(i_fnl));
v3=m(1);
a3=m(2);

%determine vel slope at inital part of record (t1 - t2)
% . . . . . . . . . . . . . . . . . . . . . . . . . . 
i_int = find(and(t>=t1, t<t2)); %times between t1 and t2
%least squares fit
% G = [ones(size(v_inp(i_int))) t(i_int)];
% m=lsqlin(G,v_inp(i_int));
% v1=m(1);
% a1=m(2);
G = [t(i_int)-t1];
m=lsqlin(G,v_inp(i_int));
v1=0;
a1=m(1);

%fit acceleration in middle  part of record (t2 - t3)
% . . . . . . . . . . . . . . . . . . . . . . . . . . 
i_mid = find(and(t>=t2, t<t3)); %times between t2 and t3
%baseline velocities at t2 and t3
v_t2 = v1 + a1*(t2-t1); %based on t1-t2 segment
v_t3 = v3 + a3*t3; %based on t2-t3 segment
%acc in middle part
a2 = (v_t3-v_t2)/(t3-t2);

% Time-series correction
%---   ---   ---   ---   ---   ---   ---   ---   ---
%remove acc baseline acc
a        = a_inp;
a(i_int) = a(i_int)-a1;
a(i_mid) = a(i_mid)-a2;
a(i_fnl) = a(i_fnl)-a3;
%compute vel and disp
v=cumtrapz(t,a);
d=cumtrapz(t,v);

%correction
a_corr = a_inp-a; 
v_corr = cumtrapz(t,a_corr);
d_corr = cumtrapz(t,v_corr);

end