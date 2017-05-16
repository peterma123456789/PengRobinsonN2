function [ p ] = getPfromTandRho( T,rho )
% getPfromTandRho Compute pressure given temperature and density

[a,b,R,dadT,d2adT2] = getThermo(T);

v = 1./rho;
p = R*T./(v-b) - a./(v.^2+2*v*b-b^2);

end

