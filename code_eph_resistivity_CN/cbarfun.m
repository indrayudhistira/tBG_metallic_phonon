function Ans = cbarfun(teV)

if nargin < 1 || isempty(teV), teV = 3.4; end

hbar = 1.0545718e-34;
e = 1.60217662e-19;

c = 1.623e4;
b = 2.46e-10;

t = teV * e;

Ans = (2 / sqrt(3)) * hbar * c / (b * t);

