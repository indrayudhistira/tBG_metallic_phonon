function Ans = sigmaephCN_highT_intraband_avg_iso(Texp, nexp, deg, decompose_mode, teV, betaAeV, vbar)

if nargin < 4 || isempty(decompose_mode), decompose_mode = 1; end %intravalley only
if nargin < 5 || isempty(teV), teV = 2.7; end
if nargin < 6 || isempty(betaAeV), betaAeV = 3.6; end
if nargin < 7 || isempty(vbar), vbar = vbarfun(deg); end

sigmaxx = sigmaephCN_highT_intraband_iso(Texp, nexp, deg, 1, decompose_mode, teV, betaAeV, vbar);
sigmayy = sigmaephCN_highT_intraband_iso(Texp, nexp, deg, 2, decompose_mode, teV, betaAeV, vbar);

Ans = sqrt(sigmaxx .* sigmayy);