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

for i=1:his.npt-1
    intaccSq(1)   = 0 ;
    intvelSq(1)   = 0 ;
    intdisSq(1)   = 0 ;
    intaccSq(i+1) = intaccSq(i) + (accSq(i+1)+accSq(i)) * dt/2 ;
    intvelSq(i+1) = intvelSq(i) + (velSq(i+1)+velSq(i)) * dt/2 ;
    intdisSq(i+1) = intdisSq(i) + (disSq(i+1)+disSq(i)) * dt/2 ;
end

his.Id_ofT = intdisSq ;
his.Iv_ofT = intvelSq ;
his.Ia_ofT = pi / (2*9.81) * intaccSq ;

his.Id     = max (his.Id_ofT) ;
his.Iv     = max (his.Iv_ofT) ;
his.Ia     = max (his.Ia_ofT) ;

his.IdNorm_ofT = his.Id_ofT / his.Id;
his.IvNorm_ofT = his.Iv_ofT / his.Iv;
his.IaNorm_ofT = his.Ia_ofT / his.Ia;


% calculate D5_95
his.t5  = his.time(min(find(his.IaNorm_ofT >= 0.05))) ;
his.t95 = his.time(min(find(his.IaNorm_ofT >= 0.95)));
his.t75 = his.time(min(find(his.IaNorm_ofT >= 0.75))) ;

his.D5_95 = his.t95 - his.t5 ;
his.D5_75 = his.t75 - his.t5 ;

end
