function [g1,g2] = exp_sigmoid_erf(varsigma,v0,v,v_var)

% This is the expectation of a Gaussian through the erf sigmoid
% It is also the first term in the expresion for E[xg(x)]
g1 = 0.5.*erf((v - v0) ./ sqrt(2*(varsigma + v_var))) + 0.5;

% This is the 2nd term in the expression for E[xg(x)] (cov * g2)
g2 = 1./(sqrt(2*pi*(varsigma + v_var))) .* ...
    exp(-(v - v0).^2 ./ (2*(varsigma + v_var)));