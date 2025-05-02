function Ans = invtauephtildeCN_highT_intraband_intravalley_iso(r, phi0, deg, teV, betaAeV, vbar)

%Units
%E, r = EVHS
%deg = degrees
%phi0 = radian
%Ans = hbar / (kB T)
%betaAeV = eV
%betaA, t = Joule
%ktheta = 2 pi / b
%b = m
%c = m / s
%mus = kg / m^2

if nargin < 4 || isempty(teV), teV = 2.7; end
if nargin < 5 || isempty(betaAeV), betaAeV = 3.6; end
if nargin < 6 || isempty(vbar), vbar = vbarfun(deg); end

%hbar = 1.0545718e-34;
e = 1.60217662e-19;

rad = (deg / 180) * pi;
gamma = 1 ./ (2 * tan(rad / 2));

betaA = betaAeV * e * vbar * gamma;
c = 1.623e4;
mus = 7.66e-7;
b = 2.46e-10;
t = teV * e;

ktheta = (4 / 3) * sin(rad / 2);

J = @(phi) (1 / 16) * r ./ sqrt((1 + r .* cos(phi)) .^ 2 + (r .* sin(phi)) .^ 2);
%J = @(phi) r; %DIRAC

F = @(phi) (1 + cos(phi)) / 2;

transportfactor = @(phi) 1 - cos(phi); %WATCH !

integrand = @(phi) J(phi) .* F(phi - phi0) .* transportfactor(phi - phi0);
angular = integral(integrand, -pi, pi, 'ArrayValued', true, 'RelTol', 1e-3);
Ans = (64 / sqrt(3)) * ((betaA ^ 2) / (mus * t * c ^ 2 * b ^ 2)) * (ktheta / vbar) * angular;