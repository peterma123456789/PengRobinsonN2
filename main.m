%%

clear all, close all, clc;

T = linspace(100,300,1000);
p = 40e5*ones(size(T));

rho = getRhofromTandP(T,p);
cp = getCpfromTandRho(T,rho);
sos = getSosfromTandRho(T,rho);

figure,plot(T,sos)







