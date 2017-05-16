function [ cp ] = getCpfromTandRho( T,rho )
% getSos Compute speed of sound given primitive variables. NASA polynomial
% for N2 is used.

% nasa polynomial for N2 (300-1000K)
an = [3.531005280E+00,-1.236609870E-04,-5.029994370E-07,2.435306120E-09,...
     -1.408812350E-12,-1.046976280E+03,2.967474680E+00];

v = 1./rho;

[a,b,R,dadT,d2adT2] = getThermo(T);
cp_ideal = R*(an(1) + an(2)*T + an(3)*T.^2 + an(4)*T.^3 + an(5)*T.^4);
cv_ideal = cp_ideal - R;

dpdT = R./(v-b) - dadT./(v.^2+2*v.*b-b.^2);
dpdv = -R.*T./(v-b).^2.*(1-2*a.*((R.*T.*(v+b).*((v.^2+2*v.*b-b.^2)./(v.^2-b.^2)).^2).^(-1)));
K1 = 1/sqrt(8)./b.*log((v+(1-sqrt(2)).*b)./(v+(1+sqrt(2)).*b));
cp = cv_ideal - T.*(dpdT).^2./dpdv - K1.*T.*d2adT2;

end

