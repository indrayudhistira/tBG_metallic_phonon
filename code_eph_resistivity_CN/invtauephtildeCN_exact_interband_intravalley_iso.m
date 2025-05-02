function Ans = invtauephtildeCN_exact_interband_intravalley_iso(r, phi0, mu, T, deg, teV, betaAeV, vbar)

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

%rprime = @(phi) r; %QE
lambda = -1;
s = 1; %absorption
rprime_helper = @(phi) rprimehybridfun(r, phi0, z, phi, s, lambda);
rprime = memoize(rprime_helper);

J = @(phi) (1 / 16) * rprime(phi) ./ sqrt((1 + rprime(phi) .* cos(phi)) .^ 2 + (rprime(phi) .* sin(phi)) .^ 2);
%J = @(phi) rprime(phi) / 16; %DIRAC

Bprime = @(phi) (1 / 8) * (1 + rprime(phi) .* cos(phi));
Cprime = @(phi) (1 / 8) * rprime(phi) .* sin(phi);
Aprime = @(phi) sqrt(Bprime(phi) .^ 2 + Cprime(phi) .^ 2);
kxprimesq = @(phi) Aprime(phi) - Bprime(phi);
kyprimesq = @(phi) Aprime(phi) + Bprime(phi);

%CNfactor = @(l, phi) 1 ./ (8 * Aprime(l, phi));
%qr is q times r
qr = @(phi) 4 .* sqrt(kxprimesq(phi) + kxsq + kyprimesq(phi) + kysq ...
                         - 2 * sign(sin(phi0) .* sin(phi)) .* sqrt(kxprimesq(phi) .* kxsq) ...
                         - 2 * sqrt(kyprimesq(phi) .* kysq));   

bosehighT = @(phi) T ./ (z * qr(phi));
bose = @(phi) 1 ./ (exp(1 ./ bosehighT(phi)) - 1);
fermi = @(phi) 1 ./ (exp(rprime(phi) / T - mu) + 1);

occabsorption = @(phi) fermi(phi) + bose(phi);

    function Ans1 = occfactor(phi)        
        Ans1 = occabsorption(phi) ./ bosehighT(phi); %ratio
        Ans1(Ans1 == Inf | isnan(Ans1)) = 1; %for r approx 0 OR r = 0
    end

integrand = @(phi) (1 / 2) * J(phi) .* occfactor(phi) .* F(phi - phi0) .* transportfactor(phi - phi0);

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