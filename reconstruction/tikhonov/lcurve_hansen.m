% Based on the regtools by:
%
% Per Christian Hansen, DTU Compute, April 14, 2003.

% Reference: A. N. Tikhonov & V. Y. Arsenin, "Solutions of Ill-Posed
% Problems", Wiley, 1977.

function [rho, eta, lambdas] = lcurve_hansen(U, sm, b, logLambdaLimits, numSamples)

[m,n] = size(U);
[p,ps] = size(sm);

beta = U'*b;
beta2 = norm(b)^2 - norm(beta)^2;

if (ps==1)
  s = sm; beta = beta(1:p);
else
  s = sm(p:-1:1,1)./sm(p:-1:1,2); beta = beta(p:-1:1);
end

xi = beta(1:p)./s;
xi(isinf(xi)) = 0;

lambdas = logspace(logLambdaLimits(1), logLambdaLimits(2), numSamples)';
eta = NaN(numSamples,1);
rho = eta;
s2 = s.^2;

for i = 1:numSamples
    f = s2./(s2 + lambdas(i)^2);
    eta(i) = norm(f.*xi);
    rho(i) = norm((1-f).*beta(1:p));
end

if m > n && beta2 > 0
	rho = sqrt(rho.^2 + beta2);
end

end