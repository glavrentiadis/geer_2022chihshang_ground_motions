function [vel_pulse_disp] = VelPulseMP03Disp(t, A,f_p,nu,gamma,t0)
%VelPulseMP03Disp displacement time history of MP03 velocity pulse
    
%velocity pulse
vel_pulse      = VelPulseMP03(t, A,f_p,nu,gamma,t0);
%permanent slip
vel_pulse_disp = cumtrapz(t, vel_pulse);

end