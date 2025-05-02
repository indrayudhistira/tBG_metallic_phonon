function output=Fmu(x) 

% This code allows you to convert from charge density n to a chemical
% potential mu(T). 

% x= T/TF. To get mu, simply multiply Fmu(x) by
% EF=sqrt(hbar*vf*sqrt(pi*abs(n)))*sign(n). Note that TF = abs(EF)/kBT. It
% is never negative.re

x = abs(x);

g = @(x) ( 1+ erf(10.*(x-0.5)) ) ./2 ; 

gbar = @(x) erfc(10.*(x-0.5))./2 ;

output = gbar(x).*(1-pi.^2.*x.^2./6) + g(x)./(4.*log(2).*x); 

 output( x  < 1e-8  ) = 1; % For x < 1e-9 and below, Fmu starts to diverge
% in the positive direction. But be mindful that we are doing this step
% approximation which will wash out all detail at x < 1e-9.

%  output( x  == 0  ) = 1; % This is what Indra does in his code.

output( x == Inf) = 0;

end