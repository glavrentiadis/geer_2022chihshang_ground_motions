function [Sa,varargout] = rspec(t,acc,T,damp)
%
% Usage: Sa = rspec(t,acc,T,damp)
%        [Sa,Sv] = rspec(t,acc,T,damp)
%        [Sa,Sv,Sd] = rspec(t,acc,T,damp)
%
% where: Sa = absolute spectral acceleration
%        Sv = relative spectral velocity
%        Sd = relative spectral displacement
%
%        t = vector of time values
%        acc = vector of ground surface acceleration
%        T = vector of natural periods
%        damp = decimal value of spectral damping
%
% Script to compute response spectra from acceleration time histories.
% This script is a Matlab "translation" of subroutine rdrvaa.for by:
%
% Boore,D. (2008), "TSPP---A Collection of FORTRAN Programs for 
% Processing and Manipulating Time Series," U.S. GEOLOGICAL SURVEY 
% OPEN-FILE REPORT 08--xxx, Version 1.0, February 2008
%
% This function is based on the method described in:
%
% Nigam, N.C., and P.C. Jennings. (1969). "Calculation of Response Spectra
% from Strong-Motion Earthquake Records." Bulletin of the Seismological
% Society of America. 59(2), 909-922.

%Reset timing to zero start
t = t - t(1);

%Length of time series
na = length(acc);

%Circular frequency
omega = 2*pi./T;

%Initial values of displacement and velocity
d0 = 0;
v0 = 0;

%Initialize vectors for spectral quantities
Sd = zeros(length(T),1);
Sv = zeros(length(T),1);
Sa = zeros(length(T),1);

%Constants
dt = t(2)-t(1);
d2 = sqrt(1-damp^2);
n1 = na - 1;

for j = 1:length(T)
    
    omt = omega(j)*dt;
    omd = omega(j)*d2;    
    omdt = omd*dt;
    om2 = omega(j)^2;
    
    bom = damp*omega(j);
    bomt = damp*omt;    
    d3 = 2*bom; % for Sa
    
    c1 = 1/om2;
    c2 = 2*damp/(om2*omt);
    c3 = c1 + c2;
    c4 = 1/(omega(j)*omt);
    
    ss = sin(omdt);
    cc = cos(omdt);   
    
    ee = exp(-bomt);
    ss = ss*ee;
    cc = cc*ee;
    
    s1 = ss/omd;
    s2 = s1*bom;
    s3 = s2 + cc;
    s4 = c4*(1-s3);
    s5 = s1*c4 + c2;   
    
    a11 = s3;
    a12 = s1;
    a21 = -om2*s1;
    a22 = cc - s2;
    
    b11 = s3*c3 - s5;
    b12 = -c2*s3 + s5 - c1;
    b21 = -s1 + s4;
    b22 = -s4;
    
    rd = 0;
    rv = 0;
    aa = 0;
    
    y = -d0;
    ydot = -v0; % minus because forcing function is -acc

    for i = 1:n1
    
        y1 = a11*y + a12*ydot + b11*acc(i) + b12*acc(i+1);   
        ydot = a21*y + a22*ydot + b21*acc(i) + b22*acc(i+1);
    
        y = y1 ;   % y is the oscillator output at time corresponding to index i

        z = abs(y);
        if z > rd
            rd = z;
        end
        
        z1 = abs(ydot) ; % for rv
        if z1 > rv
            rv = z1; % for rv
        end
        
        ra = -d3*ydot -om2*y1; % for aa
        
        z2 = abs(ra); % for aa
        if z2 > aa
            aa = z2; % for aa
        end
    end
    
    Sd(j) = rd;
    Sv(j) = rv;
    Sa(j) = aa;
    
end

if (nargout >= 2)
    varargout{1} = Sv;
end

if (nargout == 3)
    varargout{2} = Sd;
end
