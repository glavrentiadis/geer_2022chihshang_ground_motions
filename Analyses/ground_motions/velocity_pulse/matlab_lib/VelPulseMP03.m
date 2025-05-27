function [vel_pulse] = VelPulseMP03(t, A,f_p,nu,gamma,t0)
%VelPulseMP03 velocity pulse functional form from:
% G.P. Mavroeidis and A.S. Papageorgiou (2003), A Mathematical Representation of
% Near-Fault Ground Motions, Bulletin of the Seismological Society of
% America.

%vel pulse
vel_pulse = A/2*(1 + cos(2*pi*f_p/gamma * (t-t0))) .* cos(2*pi*f_p*(t-t0) + nu);   %unmasked vel pulse
t_mask    = and( t >= t0-gamma/(2*f_p), t <= t0+gamma/(2*f_p) );                         %mask
vel_pulse = vel_pulse .* t_mask;                                                        %masked velocity pulse

end
