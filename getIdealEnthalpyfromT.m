function [ h ] = getIdealEnthalpyfromT( T )
%GETIDEALENTHALPYFROMT Summary of this function goes here
%   Detailed explanation goes here

% nasa polynomial for N2 (300-1000K)
an = [3.531005280E+00,-1.236609870E-04,-5.029994370E-07,2.435306120E-09,...
     -1.408812350E-12,-1.046976280E+03,2.967474680E+00];

[a,b,R,dadT,d2adT2] = getThermo(T);

h = R*T.*(an(1) + an(2)*T/2 + an(3)*T.^2/3 + an(4)*T.^3/4 + an(5)*T.^4/5 + an(6)./T);

end

