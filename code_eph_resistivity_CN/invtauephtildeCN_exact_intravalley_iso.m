function Ans = invtauephtildeCN_exact_intravalley_iso(r, phi0, mu, T, deg, teV, betaAeV, vbar)

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

if vbar >= cbarfun(teV)
    Ans = invtauephtildeCN_exact_intraband_intravalley_iso(r, phi0, mu, T, deg, teV, betaAeV, vbar);
else    
    Ans = invtauephtildeCN_exact_interband_intravalley_iso(r, phi0, mu, T, deg, teV, betaAeV, vbar);
end 