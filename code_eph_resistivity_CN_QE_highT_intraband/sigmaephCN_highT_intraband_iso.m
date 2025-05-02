function Ans = sigmaephCN_highT_intraband_iso(Texp, nexp, deg, direction_mode, decompose_mode, teV, betaAeV, vbar)
    if nargin < 5 || isempty(decompose_mode), decompose_mode = 1; end %intravalley only
    if nargin < 6 || isempty(teV), teV = 2.7; end
    if nargin < 7 || isempty(betaAeV), betaAeV = 3.6; end
    if nargin < 8 || isempty(vbar), vbar = vbarfun(deg); end
    
    onesmatrix = ones(size(Texp .* nexp .* deg .* vbar));
    Ans = arrayfun(@(TTexp, nnexp, ddeg, vvbar) sigmaephCN_highT_intraband_iso_scalar(TTexp, nnexp, ddeg, direction_mode, decompose_mode, teV, betaAeV, vvbar), onesmatrix .* Texp, onesmatrix .* nexp, onesmatrix .* deg, onesmatrix .* vbar);
end

function Ans = sigmaephCN_highT_intraband_iso_scalar(Texp, nexp, deg, direction_mode, decompose_mode, teV, betaAeV, vbar)

%Castro Neto model for twisted bilayer graphene

%Units
%Texp, TVHS = Kelvin
%nexp, nVHS = 10^10 cm^-2
%deg = degrees
%mu = kB T
%EVHSeV = eV
%E = EVHS
%T = TVHS
%n = nVHS
%kx, ky = ktheta
%vx = vF / 4

%mode: 1 for sigmaxx; 2 for sigmayy

tic;

kB = 1.38064852e-23;
e = 1.60217662e-19;

EVHSeV = EVHSfunCN(deg, teV, vbar);
TVHS = EVHSeV * e / kB;

T = Texp / TVHS;
n = nexp / nVHSfunCN(deg);

%Only dimensionless T and n enter calculation from now onwards

mu = mutildefunCN(T, n);
if isnan(mu)
    Ans = NaN;
    return;
end

y = @(r, phi) sqrt((1 + r .* cos(phi)) .^ 2 + (r .* sin(phi)) .^ 2);
J = @(r, phi) (1 / 16) * r ./ y(r, phi);

kxsq = @(r, phi) (-(1 + r .* cos(phi)) + y(r, phi)) / 8;
kysq = @(r, phi) (1 + r .* cos(phi) + y(r, phi)) / 8;

vxsq = @(r, phi) 64 * kxsq(r, phi) .* ((y(r, phi) + 1) ./ r) .^ 2;
vysq = @(r, phi) 64 * kysq(r, phi) .* ((y(r, phi) - 1) ./ r) .^ 2;

% E = @(kx, ky, s) s * sqrt(kx .^ 2 + ky .^ 2); %DIRAC
% vxsq = @(r, phi) 1; %DIRAC
% vysq = @(r, phi) 1; %DIRAC

switch direction_mode
    case 1
        visq = vxsq;
    case 2
        visq = vysq;
    otherwise
        error('Wrong selection');
end

invtaueph_intravalley = @(r, phi) invtauephtildeCN_highT_intraband_intravalley_iso(r, phi, deg, teV, betaAeV, vbar);
invtaueph_intervalley = @(r, phi) invtauephtildeCN_highT_intraband_intervalley_iso(r, phi, deg, teV, betaAeV, vbar);

switch decompose_mode
    case 1
        invtaueph = invtaueph_intravalley;
    case 2
        invtaueph = invtaueph_intervalley;
    case 3
        invtaueph = @(r, phi) invtaueph_intravalley(r, phi) + invtaueph_intervalley(r, phi);
    otherwise    
        error('Wrong selection');        
end
taueph = @(r, phi) 1 ./ invtaueph(r, phi);

%s = +1 for e; -1 for h
integrandeh = @(r, phi, s) J(r, phi) .* visq(r, phi) .* taueph(r, phi) ./ (4 * (cosh((s .* r / T - mu) / 2)) .^ 2);
%integrand = @(x, phi) J(abs(x) * T, phi) .* visq(abs(x) * T, phi) .* taueph(abs(x) * T, phi) ./ (4 * (cosh((x - mu) / 2)) .^ 2);

g = 8;

lim = 20;
Ans = (1 / pi) * (g / 2) * (1 / T) * integral2(@(x, phi) integrandeh(x * T, phi, 1) + integrandeh(x * T, phi, -1), max(0, abs(mu) - lim), abs(mu) + lim, -pi, pi, 'RelTol', 1e-4);
%Ans = (1 / pi) * (g / 2) * (1 / T) * integral2(@(x, phi) integrandeh(abs(x) * T, phi, sign(x)), mu - lim, mu + lim, -pi, pi, 'RelTol', 1e-4);

toc;

end