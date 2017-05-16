function [ a,b,R,dadT,d2adT2 ] = getThermo( T )
% getThermo Compute necessary thermodyamic parameters for the cubic
% equation of state (N2 property is assumed).

% N2
MW = 28.0134e-3;
Tc = 126.19;
pc  = 3.3958e+6;
rhoc = 313.3;
omega  = 0.03720;

c = 0.37464 + 1.54226*omega - 0.26992*omega^2;

R = 8.314/MW;
a = 0.457236*(R*Tc)^2/pc*(1+c*(1-sqrt(T/Tc))).^2;
b = 0.077796*R*Tc./pc;
G = c*sqrt(T/Tc)./(1+c*(1-sqrt(T/Tc)));
dadT = -1./T.*a.*G;
d2adT2 = 0.457236*R^2./T/2*c*(1+c)*Tc/pc.*sqrt(Tc./T);

end

