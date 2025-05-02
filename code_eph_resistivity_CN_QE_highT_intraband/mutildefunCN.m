function Ans = mutildefunCN(TbyTVHS, nbynVHS)
    onesmatrix = ones(size(TbyTVHS .* nbynVHS));
    Ans = arrayfun(@(TTbyTVHS, nnbynVHS) mutildefunCNscalar(TTbyTVHS, nnbynVHS), onesmatrix .* TbyTVHS, onesmatrix .* nbynVHS);
end

function Ans = mutildefunCNscalar(T, n)

%Do not use ellipticK, somehow it is not accurate for big T
%Use integral2 instead

%Units
%T, TFDirac = TVHS
%n, ne_min_nh = nVHS
%Ec, E = EVHS
%kx, ky = ktheta
%ktheta = 2 * pi / a
%TstarDirac = TFDirac (T in units of TFDirac)
%mutilde = kB T
%Ans is mu / (kB T)

if T == 0
    warning('Use zero temperature code instead !');
    Ans = n / T;    
    return;
end
if abs(n) > 1, warning('CN:largeDensity', 'n is larger than nVHS. Castro Neto model can''t be trusted !'); end

fermi = @(x) 1 ./ (exp(x) + 1);

%dos = @(r) (r / 2) .* (ellipticK(4 * r ./ (1 + r) .^ 2) ./ (1 + r) + ellipticK(-4 * r ./ (1 - r) .^ 2) ./ (1 - r));
%integrand = @(r, mutilde, s) dos(r) .* fermi(r / T - s * mutilde);

y = @(r, phi) sqrt((1 - r .* cos(phi)) .^ 2 + (r .* sin(phi)) .^ 2);
J = @(r, phi) (1 / 16) * r ./ y(r, phi);
integrand2 = @(r, phi, mutilde, s) 4 * J(r, phi) .* fermi(r / T - s * mutilde);

%+1 for ne; -1 for nh
lim = 20;
%ne_min_nh = @(mu) integral(@(r) integrand(r, mu, 1) - integrand(r, mu, -1), 0, (abs(mu) + lim) * T, 'RelTol', 1e-3);
%ne_min_nh = @(mu) integral(@(r) integrand(r, mu, 1) - integrand(r, mu, -1), 0, Inf, 'RelTol', 1e-3);

%ne_min_nh = @(mu) integral2(@(r, phi) integrand2(r, phi, mu, 1) - integrand2(r, phi, mu, -1), 0, (abs(mu) + lim) * T, -pi, pi, 'RelTol', 1e-3);
%ne_min_nh = @(mu) integral2(@(r, phi) integrand2(r, phi, mu, 1) - integrand2(r, phi, mu, -1), 0, Inf, -pi, pi, 'RelTol', 1e-3);

ne_min_nh = @(mu) T * integral2(@(x, phi) integrand2(x * T, phi, mu, 1) - integrand2(x * T, phi, mu, -1), 0, abs(mu) + lim, -pi, pi, 'RelTol', 1e-3);
%ne_min_nh = @(mu) T * integral2(@(x, phi) integrand2(x * T, phi, mu, 1) - integrand2(x * T, phi, mu, -1), 0, Inf, -pi, pi, 'RelTol', 1e-3);

fun = @(mu) arrayfun(@(mmu) ne_min_nh(mmu) - n, mu);

%guess for mu
TFDirac = 2 * sign(n) * sqrt(abs(n) / pi);
TstarDirac = T / TFDirac;
mutildeDirac = Fmu(TstarDirac) / TstarDirac;

%Ans = mutildeDirac; %DIRAC

options = optimoptions('fsolve', 'Display', 'none');
[temp, ~, flag] = fsolve(fun, mutildeDirac, options); %faster than fzero

if flag > 0 %fsolve compute succesfully
    Ans = temp;
else
    warning('fsolve fail to obtain accurate solution, return NaN instead.')
    Ans = NaN;
end

end