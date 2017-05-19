function [ h ] = getEnthalpyfromTandRho( T,rho )
%GETENERGYFROMTANDRHO Summary of this function goes here
%   Detailed explanation goes here

p = getPfromTandRho(T,rho);
e = getEnergyfromTandRho(T,rho);

h = e + p./rho;

end

