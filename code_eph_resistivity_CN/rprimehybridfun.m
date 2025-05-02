function Ans = rprimehybridfun(r, phi0, z, phi, s, lambda)

B = @(r, phi0) (1 / 8) * (1 + r .* cos(phi0));
C = @(r, phi0) (1 / 8) * r .* sin(phi0);
A = @(r, phi0) sqrt(B(r, phi0) .^ 2 + C(r, phi0) .^ 2);
kxsq = @(r, phi0) A(r, phi0) - B(r, phi0);
kysq = @(r, phi0) A(r, phi0) + B(r, phi0);

Bprime = @(x) (1 / 8) * (1 + x .* cos(phi));
Cprime = @(x) (1 / 8) * x .* sin(phi);
Aprime = @(x) sqrt(Bprime(x) .^ 2 + Cprime(x) .^ 2);
kxprimesq = @(x) Aprime(x) - Bprime(x);
kyprimesq = @(x) Aprime(x) + Bprime(x);

%qr is q times r or q / (k_theta / 4)
qr = @(r, phi0, x) 4 .* sqrt(kxprimesq(x) + kxsq(r, phi0) + kyprimesq(x) + kysq(r, phi0) ...
                         - 2 * sign(sin(phi0) .* sin(phi)) .* sqrt(kxprimesq(x) .* kxsq(r, phi0)) ...
                         - 2 * sqrt(kyprimesq(x) .* kysq(r, phi0)));   

%rhs = @(x) lambda * r + s * z * qr(x);
%rhs = @(x) lambda * r + s * z * abs(qr(x));
rhs = @(r, phi0, x) abs(lambda * r + s * z * abs(qr(r, phi0, x)));

%Self-consistent solution

%Initial guess
x = r;
x1 = x;
x2 = 0.1 + 1.4 * r;

maxiter = 10; %10 24 96
j = 0;
reltol = 1e-2; %1e-2
abstol = 1e-2;
epsilon = .05;
while true
    j = j + 1;
    %index = abs((x - x2) ./ (x + 1e3)) > reltol | abs(x - x2) > abstol;
    index = abs(x - x2) > max(abstol, reltol * abs(x2));
    
    if ~any(any(index)), break; end
    
    x1(index) = rhs(r(index), phi0(index), x(index));    
    x2(index) = rhs(r(index), phi0(index), x1(index));
    
    denom = (x2 - x1) - (x1 - x);
    
    index2 = abs(denom) > epsilon;
        
    x(index & index2) = abs(x2(index & index2) - ((x2(index & index2) - x1(index & index2)) .^ 2) ./ denom(index & index2));
    x(index & ~index2) = x2(index & ~index2);
    x2(index & ~index2) = x1(index & ~index2);
    
    %Break if too long, likely no solution
    if j == maxiter
        %warning('maximum iteration reached');
%         index = abs((x - x2) ./ (x + 1e-3)) > reltol | abs(x - x2) > abstol;
%         disp(['aitken: maximum iteration reached for ' sprintf('%.2f', 100 * sum(sum(index == true)) / numel(x)) ' % of all cases']);
        break;
    end
end

maxiter = 50; %50
j = 0;
while true
	j = j + 1;
    %index = abs((x - x2) ./ (x + 1e-3)) > reltol | abs(x - x2) > abstol;
    index = abs(x - x2) > max(abstol, reltol * abs(x2));
    
    if ~any(any(index)), break; end
    
    x2(index) = x(index);
    %x(index) = rhs(r(index), phi0(index), x(index));
    x(index) = (x(index) + rhs(r(index), phi0(index), x(index))) / 2; %converge faster
    
    if j == maxiter
        %warning('maximum iteration reached');
%         index = abs((x - x2) ./ (x + 1e-3)) > reltol | abs(x - x2) > abstol;
%         disp(['maximum iteration reached for ' sprintf('%.2f', 100 * sum(sum(index == true)) / numel(x)) ' % of all cases']);
        break;
    end
end

fun = @(r, phi0, x) x - rhs(r, phi0, x);
%     function Ans1 = fun(r, phi0, x)
%         numel(r)
%         Ans1 = x - rhs(r, phi0, x);
%     end

options = optimset('TolX', 1e-2);
% [x, ~, flag] = fzero(fun, x0, options);
x(index) = arrayfun(@(rr, pphi0) fzero(@(x) fun(rr, pphi0, x), rr, options), r(index), phi0(index));

%     function [x, flag] = solver(r, phi0)
%         x = zeros(size(r));
%         flag = zeros(size(r));
%         
%         for j = 1:numel(r)
%             [x(j), ~, flag(j)] = fzero(@(x) fun(r(j), phi0(j), x), r(j), options);          
%         end    
%     end    
% 
% [x(index), flag] = solver(r(index), phi0(index));
% 
% sum(sum(flag < 1))

Ans = x;

end