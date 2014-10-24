function g = mean_sigmoid_JR(varsigma,v0,v,v_var)

g = 0.5*erf((v - v0) / sqrt(2*(varsigma^2 + (1/100)^2*v_var))) + 0.5;
