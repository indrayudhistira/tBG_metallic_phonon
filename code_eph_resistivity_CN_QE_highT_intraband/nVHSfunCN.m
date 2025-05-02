function Ans = nVHSfunCN(deg)

%Unit
%bAngstrom = Angstrom
%ktheta = 2 * pi / b
%Ans = 10^10 cm^-2

g = 8;
bAngstrom = 2.46; %graphene lattice constant

rad = (deg / 180) * pi;
ktheta = (4 / 3) * sin(rad / 2);

%The (1 / 2) is from evaluation of the dimensionless integral numerically
Ans = ((1e3 / bAngstrom) ^ 2) * (g / 2) * (ktheta .^ 2) * (1 / 2);