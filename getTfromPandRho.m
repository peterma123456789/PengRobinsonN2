function [ T ] = getTfromPandRho( p,rho )
% getTfromPandRho Compute temperature given pressure and density.

CRIT = 1.0e-6;

v = 1./rho;
T = 300*ones(1,length(rho)); % initial guess
p_n = getPfromTandRho(T,rho);
diff = p_n - p;

while (max(abs(diff)) > CRIT)
    [a,b,R,dadT,d2adT2] = getThermo(T);
    dpdT = R./(v-b) - 1./(v.^2+2*v*b-b^2).*dadT;
    T = T - diff./dpdT;
    p_n = getPfromTandRho(T,rho);
    diff = p_n - p;
end

end

