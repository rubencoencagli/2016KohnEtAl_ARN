function a = getLognAmplitudes(neurons, kappa)
% Sample a set of tuning curve amplitudes from a lognormal distribution.
%   a = getLognAmplitudes(neurons, kappa) returns amplitudes for given
%   number of neurons with amoplitude variability kappa = Var[sqrt(a)].
%
% AE 2011-01-28

[m, s] = getLognParams(kappa);
a = exp(m + randn(neurons, 1) * s);


function [m, s] = getLognParams(kappa)
% Parameters of a lognormal distribution with 
%   * E[x] = 1 
%   * Var[sqrt(x)] = kappa

v2 = 0.1;
while getKappa(v2) < kappa
    v2 = 2 * v2;
end
v1 = 0;
while v2 - v1 > 1e-4
    v = (v1 + v2) / 2;
    k = getKappa(v);
    if k > kappa
        v2 = v;
    else
        v1 = v;
    end
end
[m, s] = getMS(v);


function kappa = getKappa(v)
% kappa = Var[sqrt(x)] for a lognormal with E[x] = 1 and Var[x] = v

[m, s] = getMS(v);
b = 0.01;
x = 0:b:100;
p = mylognpdf(x, m, s);
kappa = sum(p .* x) * b - (sum(p .* sqrt(x)) * b)^2;


function [m, s] = getMS(v)
% Parameters of a lognormal with E[x] = 1 and Var[x] = v

m = log(1) - 1/2 * log(v + 1);
s = sqrt(log(v + 1));


function p = mylognpdf(x, m, s)
% Lognormal probability density function (to circumvent stats toolbox)

p = exp(-(log(x) - m).^2 / 2 / s^2) ./ x / sqrt(2*pi*s^2);
p(x < eps) = 0;
