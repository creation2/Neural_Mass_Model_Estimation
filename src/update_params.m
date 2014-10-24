function params = update_params(est_params,N_p,N_syn,params)
% Updates the parameters used by column_model_JR to the estimated
% parameters
%*****INPUTS********
%parameters - estimated values of the parameters
%N_p - numbers of parameters to estimate
%N_syn - number of synapses in the model

%*******************************
%params.constant_input(1) = est_params(end);

alpha = params.alpha;
alpha(1:7:end) = est_params(1:4);
alpha(2:7:end) = est_params(21:24);
alpha(3:7:end) = est_params(25:28);
alpha(4:7:end) = est_params(5:8);
alpha(5:7:end) = est_params(9:12);
alpha(6:7:end) = est_params(13:16);
alpha(7:7:end) = est_params(17:20);
params.alpha = alpha;

end