function [ lambda ] = getConductivityfromTandRho( T,rho )
%GETENERGYFROMTANDRHO Summary of this function goes here
%   Detailed explanation goes here

an = [3.531005280E+00,-1.236609870E-04,-5.029994370E-07,2.435306120E-09,...
     -1.408812350E-12,-1.046976280E+03,2.967474680E+00];

% N2
MW = 28.0134;
Tc = 126.19;
% pc  = 3.3958e+6;
rhoc = 313.3;
Vc = MW / rhoc * 1000;
omega  = 0.03720;
dipole = 0.0;

R = 8.314/(MW/1000);

mu_R = 131.3 * dipole / sqrt(Vc * Tc);
kappa = 0;

TStar = 1.2593*T/Tc;
Neufeld = 1.16145 * TStar.^(-0.14874) + 0.52487 * exp(-0.7732 * TStar)...
            + 2.16178 * exp(-2.43787 * TStar);
Fc = 1 - 0.2756*omega+0.059035*mu_R.^4 + kappa;

% low-pressure value
mu = 40.785 * Fc .* sqrt(MW * T) ./ (Neufeld * Vc^(2.0 / 3.0)) / 10.0e6;

% coeffs
a2 = [2.4166, -0.50924, 6.6107, 14.543, 0.79274, -5.8634, 91.089];
b2 = [0.74824, -1.5094, 5.6207, -8.9139, 0.82019, 12.801, 128.11];
c2 = [-0.91858, -49.991, 64.760, -5.6379, -0.69369, 9.5893, -54.217];
d2 = [121.72, 69.983, 27.039, 74.344, 6.3173, 65.529, 523.81];
for k = 1:length(a2)
    B_coef(k) = a2(k) + b2(k) * omega + c2(k) * mu_R^4.0 + d2(k) * kappa;
end

y = rho ./ (1000 * MW) * Vc / 6.0;
G1 = (1.0 - 0.5 * y) ./ (1.0 - y).^3;
G2 = ((B_coef(1) ./ y) .* (1.0 - exp(-B_coef(4) * y))...
        + B_coef(2) * G1 .* exp(B_coef(5) * y) + B_coef(3) * G1)...
        / (B_coef(1) * B_coef(4) + B_coef(2) + B_coef(3));

cp_ig = MW*R*(an(1) + an(2)*T + an(3)*T.^2 + an(4)*T.^3 + an(5)*T.^4);
cv_ig = cp_ig - MW*R;

alpha = cv_ig / (MW*R) - 3.0 / 2.0;
beta = 0.7862 - 0.7109 * omega + 1.3168 * omega.^2;
Zr = 2.0 + 10.5 * (T / Tc).^2;
q = 0.003586 * sqrt(Tc / (MW / 1000.0)) / Vc^(2.0 / 3.0);
psi = 1.0 + alpha * ((0.215 + 0.28288 * alpha - 1.061 * beta + 0.26665 * Zr)...
            / (0.6366 + beta * Zr + 1.061 * alpha * beta));

lambda = 31.2 * mu .* psi / (MW / 1000) .* (1.0 ./ G2 + B_coef(6) * y)...
            + q * B_coef(7) .* y.^2 .* sqrt(T / Tc) .* G2;

end

