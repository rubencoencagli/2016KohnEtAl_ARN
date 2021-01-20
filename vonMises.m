function [f, df, ddf] = vonMises(ampl, phi, theta,A,B,C)
% Von-Mises tuning function and derivatives.
%   [f, df, ddf] = vonMises(ampl, phi, theta) evaluates a von-Mises tuning
%   function and its first and second derivative at stimulus value theta.
%   Neurons are assumed to have tuning functions of identical shape but
%   different amplitudes (ampl) and preferred orientations (phi).
%
% AE 2011-02-22

alpha = A;
gamma = B;
beta = C * exp(-gamma);

psi = theta - phi;
alpha = ampl * alpha;
g = ampl .* (beta * exp(gamma * cos(psi)));

f = alpha + g;
if nargout > 1
    df = gamma * -sin(psi) .* g;
end
if nargout > 2
    ddf = gamma * (gamma * sin(psi).^2 - cos(psi)) .* g;
end
