
% note, the maximum firing rate is lumped into synaptic gain

function g = sigmoid_JR(varsigma,v0,v)

%g = 1./(1 + exp(varsigma*(v0 - v)));% - 1./(1 + exp(varsigma*v0));      % note, the maximum firing rate is lumped into synaptic gain
g = 0.5*erf((v - v0) / (sqrt(2)*varsigma)) + 0.5;

