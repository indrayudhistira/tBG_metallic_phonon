function Ans = vbarfun2(w,theta)

%Units
%alpha0 = dimensionless
%theta = degrees
%w = eV
%v0 = 1e6 m/s

thetarad = (theta / 180) * pi;


% v0 = 0.8739;
% gamma = w ./ (22.41554 * v0 * sin(thetarad / 2)); %1 / 24.6571 = (3 / (8 * pi)) * a / (hbar * v0) %1 / 22.41554 = (3 / (8 * pi)) * a / hbar
t = 2.7;
gamma = (sqrt(3) / (4 * pi)) * (w / t) ./ sin(thetarad / 2);

Ans = (1 - 3 * gamma .^ 2) ./ (1 + 6 * gamma .^ 2); %For theta >= 1st magic angle