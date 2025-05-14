function [vel,disp] = acce2veldisp(acce,dt)
%Script to compute velocity and displacement from acceleration time
%history.

%number of points
N = max(size(acce));

%constants
dt2 = 0.5*dt;
dt6 = dt^2/6;

%numerical integration
%initalization
cumvel = 0;
cumdisp = 0;
vel  = zeros(N,1);
disp = zeros(N,1);

for i = 2:N
    a1 = acce(i-1);
    a2 = acce(i);
    v1 = vel(i-1);
    cumvel = cumvel + (a1 + a2)*dt2;
    vel(i) = cumvel;
    cumdisp = cumdisp + v1*dt + (2.0*a1 + a2)*dt6;
    disp(i) = cumdisp;
end

end