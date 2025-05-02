function Ans = sigmaephCN_exact_avg_iso(Texp, nexp, deg, decompose_mode, teV, betaAeV, vbar)

if nargin < 4 || isempty(decompose_mode), decompose_mode = 1; end %intravalley only
if nargin < 5 || isempty(teV), teV = 3.4; end
if nargin < 6 || isempty(betaAeV), betaAeV = 3.6; end
if nargin < 7 || isempty(vbar), vbar = vbarfun(deg); end

sigmaxx = sigmaephCN_exact_iso(Texp, nexp, deg, 1, decompose_mode, teV, betaAeV, vbar);
sigmayy = sigmaephCN_exact_iso(Texp, nexp, deg, 2, decompose_mode, teV, betaAeV, vbar);

Ans = sqrt(sigmaxx .* sigmayy);