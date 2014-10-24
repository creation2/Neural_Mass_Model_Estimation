function g = sigmoid_erf(varsigma,v0,v)

varsigma = repmat(varsigma,1,size(v,2));
v0 = repmat(v0,1,size(v,2));
g = 0.5.*erf((v - v0) ./ (sqrt(2).*varsigma)) + 0.5;