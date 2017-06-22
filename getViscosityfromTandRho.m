function [ mu ] = getViscosityfromTandRho( T,rho )
%GETENERGYFROMTANDRHO Summary of this function goes here
%   Detailed explanation goes here

% N2
MW = 28.0134;
Tc = 126.19;
% pc  = 3.3958e+6;
rhoc = 313.3;
Vc = MW / rhoc * 1000;
omega  = 0.03720;
dipole = 0.0;

mu_R = 131.3 * dipole / sqrt(Vc * Tc);
kappa = 0;

TStar = 1.2593*T/Tc;
Neufeld = 1.16145 * TStar.^(-0.14874) + 0.52487 * exp(-0.7732 * TStar)...
            + 2.16178 * exp(-2.43787 * TStar);
Fc = 1 - 0.2756*omega+0.059035*mu_R.^4 + kappa;

% coeffs
a1 = [6.324, 1.210e-3, 5.283, 6.623, 19.745, -1.900, 24.275, 0.7972, -0.2382, 0.06863];
b1 = [50.412, -1.154e-3, 254.209, 38.096, 7.630, -12.537, 3.450, 1.1170, 0.0677, 0.3479];
c1 = [-51.680, -6.257e-3, -168.480, -8.464, -14.354, 4.985, -11.291, 0.01235, -0.8163, 0.5926];
d1 = [1189.000, 0.03728, 3898.0, 31.42, 31.53, -18.15, 69.35, -4.117, 4.025, -0.727];
for k = 1:length(a1)
    E(k) = a1(k) + b1(k) * omega + c1(k) * mu_R.^4.0 + d1(k) * kappa;
end

y = rho ./ (1000 * MW) * Vc / 6.0;
G1 = (1.0 - 0.5 * y) ./ (1.0 - y).^3;
G2 = (E(1) ./ y .* (1.0 - exp(-E(4) * y)) + E(2) * G1 .* exp(E(5) * y)...
        + E(3) * G1) ./ (E(1) * E(4) + E(2) + E(3));

mustarstar = E(7) * y.^2 .* G2...
        .* exp(E(8) + E(9) ./ TStar + E(10) ./ TStar.^2);
mustar = TStar.^0.5 ./ Neufeld .* (Fc .* (1.0 ./ G2 + E(6) * y)) + mustarstar;
 
mu = (36.644 * mustar * sqrt(MW * Tc) / Vc.^(2.0 / 3.0)) / 10.0e6;

end

