function Ans = vbarfun(deg, wbyt)

%Units
%alpha0 = dimensionless
%theta = degrees
%w = eV
%v0 = 1e6 m/s

if nargin < 2 || isempty(wbyt), wbyt = 0.1317185 / 3.4; end

rad = (deg / 180) * pi;

%w = 0.11;
% v0 = 0.8977;
% gamma = w ./ (22.41554 * v0 * sin(rad / 2)); %1 / 24.6571 = (3 / (8 * pi)) * a / (hbar * v0) %1 / 22.41554 = (3 / (8 * pi)) * a / (hbar * v0)
%t = 2.7;

%gamma = (sqrt(3) / (4 * pi)) * (w / t) ./ sin(rad / 2);
gamma = (sqrt(3) / (4 * pi)) * wbyt ./ sin(rad / 2);

Ans = (1 - 3 * gamma .^ 2) ./ (1 + 6 * gamma .^ 2); %For theta >= 1st magic angle