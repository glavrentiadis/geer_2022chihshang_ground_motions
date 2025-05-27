function [his] = husid (his, grav_acc) 
% calculate arias intensity and associated intensity parameters
% -----------------------------------------------------------------------
% arias intensity and D5_95

%default input
if nargin < 2; grav_acc = 1.; end

%constants
gravity = 9.80665; 

%gravity correction ratio from acc units to m/sec^2
grav_ratio = gravity / grav_acc;

% calculate arias intensity

dt = his.dt;
accSq = (his.acc*grav_ratio) .^ 2 ;
velSq = his.vel .^ 2 ;
disSq = his.dis .^ 2 ;

intaccSq = zeros(his.npt,1);
intvelSq = zeros(his.npt,1);
intdisSq = zeros(his.npt,1);
for i=1:his.npt-1
    intaccSq(i+1) = intaccSq(i) + (accSq(i+1)+accSq(i)) * dt/2 ;
    intvelSq(i+1) = intvelSq(i) + (velSq(i+1)+velSq(i)) * dt/2 ;
    intdisSq(i+1) = intdisSq(i) + (disSq(i+1)+disSq(i)) * dt/2 ;
end

his.Id_his = intdisSq ;
his.Iv_his = intvelSq ;
his.Ia_his = pi / (2*9.81) * intaccSq ;

his.Id     = max (his.Id_his) ;
his.Iv     = max (his.Iv_his) ;
his.Ia     = max (his.Ia_his) ;

his.IdNorm_his = his.Id_his / his.Id;
his.IvNorm_his = his.Iv_his / his.Iv;
his.IaNorm_his = his.Ia_his / his.Ia;


% calculate D5_95
his.t5  = his.time(find(his.IaNorm_his >= 0.05, 1 )) ;
his.t95 = his.time(find(his.IaNorm_his >= 0.95, 1 ));
his.t75 = his.time(find(his.IaNorm_his >= 0.75, 1 )) ;

his.D5_95 = his.t95 - his.t5 ;
his.D5_75 = his.t75 - his.t5 ;

end
