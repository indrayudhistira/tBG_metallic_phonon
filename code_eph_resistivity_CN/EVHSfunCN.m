function Ans = EVHSfunCN(deg, teV, vbar)

%Units
%teV, Ans = eV
%ktheta = 2 pi / b

if nargin < 2 || isempty(teV), teV = 3.4; end
if nargin < 3 || isempty(vbar), vbar = vbarfun(deg); end

rad = (deg / 180) * pi;
ktheta = (4 / 3) * sin(rad / 2);

Ans = (sqrt(3) * pi / 4) * teV * vbar .* ktheta;