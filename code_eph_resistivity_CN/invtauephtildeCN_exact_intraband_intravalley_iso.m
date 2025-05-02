function Ans = invtauephtildeCN_exact_intraband_intravalley_iso(r, phi0, mu, T, deg, teV, betaAeV, vbar)

%Units
%E, r = EVHS
%phi0 = radian
%mu = kB T
%T = T_VHS
%deg = degrees
%Ans = hbar / (kB T)
%betaAeV = eV
%betaA, t = Joule
%ktheta = 2 pi / b
%b = m
%c = m / s
%mus = kg / m^2

if nargin < 6 || isempty(teV), teV = 3.4; end
if nargin < 7 || isempty(betaAeV), betaAeV = 3.6; end
if nargin < 8 || isempty(vbar), vbar = vbarfun(deg); end

% %Avoid r == 1
% r(r == 1) = 1 + 1e-5; %gives problem when r == 1 and phi0 != 0

e = 1.60217662e-19;

rad = (deg / 180) * pi;
gamma = 1 ./ (2 * tan(rad / 2));

betaA = betaAeV * e * vbar * gamma;
c = 1.623e4;
mus = 7.66e-7;
b = 2.46e-10;
t = teV * e;

z = cbarfun(teV) / vbar;

ktheta = (4 / 3) * sin(rad / 2);

%For F and transportfactor, phi here is phi_{p,p'}
F = @(phi) (1 + cos(phi)) / 2; 
transportfactor = @(phi) 1 - cos(phi); %WATCH !

B = (1 / 8) * (1 + r .* cos(phi0));
C = (1 / 8) * r .* sin(phi0);
A = sqrt(B .^ 2 + C .^ 2);
kxsq = A - B;
kysq = A + B;

%rprime = @(phi, s) r; %QE
lambda = 1;
rprime_helper = @(phi, s) rprimehybridfun(r, phi0, z, phi, s, lambda);
rprime = memoize(rprime_helper);

J = @(phi, s) (1 / 16) * rprime(phi, s) ./ sqrt((1 + rprime(phi, s) .* cos(phi)) .^ 2 + (rprime(phi, s) .* sin(phi)) .^ 2);
%J = @(phi, s) rprime(phi, s) / 16; %DIRAC

Bprime = @(phi, s) (1 / 8) * (1 + rprime(phi, s) .* cos(phi));
Cprime = @(phi, s) (1 / 8) * rprime(phi, s) .* sin(phi);
Aprime = @(phi, s) sqrt(Bprime(phi, s) .^ 2 + Cprime(phi, s) .^ 2);
kxprimesq = @(phi, s) Aprime(phi, s) - Bprime(phi, s);
kyprimesq = @(phi, s) Aprime(phi, s) + Bprime(phi, s);

%CNfactor = @(l, phi) 1 ./ (8 * Aprime(l, phi));
%qr is q times r
qr = @(phi, s) 4 .* sqrt(kxprimesq(phi, s) + kxsq + kyprimesq(phi, s) + kysq ...
                         - 2 * sign(sin(phi0) .* sin(phi)) .* sqrt(kxprimesq(phi, s) .* kxsq) ...
                         - 2 * sqrt(kyprimesq(phi, s) .* kysq));   

bosehighT = @(phi, s) T ./ (z * qr(phi, s));
bose = @(phi, s) 1 ./ (exp(1 ./ bosehighT(phi, s)) - 1);
fermi = @(phi, s) 1 ./ (exp(rprime(phi, s) / T - mu) + 1);

% occabsorption = @(phi) fermi(phi, 1) + bose(phi, 1);
% occemission = @(phi) 1 - fermi(phi, -1) + bose(phi, -1);
occ = @(phi, s) (1 - s) / 2 + s * fermi(phi, s) + bose(phi, s);

    function Ans1 = occfactor(phi, s)        
        Ans1 = occ(phi, s) ./ bosehighT(phi, s); %ratio
        Ans1(Ans1 == Inf | isnan(Ans1)) = 1; %for r approx 0 OR r = 0
    end

Jxoccfactor = @(phi, s) J(phi, s) .* occfactor(phi, s);
integrand = @(phi) (1 / 2) * (Jxoccfactor(phi, 1) + Jxoccfactor(phi, -1)) .* F(phi - phi0) .* transportfactor(phi - phi0);

dphi = 1e-2; %5e-3
%angular = integral(integrand, -pi, pi);
%angular = integral(integrand, -pi+dphi, pi-dphi, 'ArrayValued', true, 'RelTol', 1e-3);

N = 120; %160
deltaphi = 2 * (pi - dphi) / N;

angular = 0;
%Midpoint rule
phi = linspace(-pi + dphi + deltaphi / 2, pi - dphi - deltaphi / 2, N);
for j = 1:N
    angular = angular + integrand(phi(j));
end    
% % Trapezoidal rule
% phi = linspace(-pi + dphi, pi - dphi, N + 1);
% for j = 1:(N + 1)
%     if j == 1 || j == N + 1
%         angular = angular + integrand(phi(j)) / 2;
%     else
%         angular = angular + integrand(phi(j));
%     end
% end
angular = angular * deltaphi;

Ans = (64 / sqrt(3)) * ((betaA ^ 2) / (mus * t * c ^ 2 * b ^ 2)) * (ktheta / vbar) * angular;

end