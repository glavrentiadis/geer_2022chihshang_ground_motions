function [acce,disp] = vel2accedisp(vel,dt)
%Script to compute acceleration and displacement from velocity time
%history.

%Number of points
N = max(size(vel));

%acceleration and displacement
acce = [0;diff(vel)/dt];
disp = dt*cumtrapz(vel);

end