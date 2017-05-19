function [ e ] = getEnergyfromTandRho( T,rho )
%GETENERGYFROMTANDRHO Summary of this function goes here
%   Detailed explanation goes here

v = 1./rho;

[a,b,R,dadT,d2adT2] = getThermo(T);

K1 = 1/sqrt(8)/b*log((v+(1-sqrt(2))*b)./(v+(1+sqrt(2))*b));

e = getIdealEnthalpyfromT(T) - R*T + (a - T.*dadT).*K1;

end

