
function y = observation_function(x,H)


y = H*x(1:size(H,2),:);