% Models a Gaussian through a neural mass column/s
%
% ************ INPUTS *************
% x - state vector of membrane voltages and derivatives (alternating
% i.e. [v1 z1 ... vN zN]
% P - covariance of x
% Q - process noise/predicition uncertainty for x
% params - struct containing model parameters
%
%************* OUTPUTS ***************
% x_out - mean of the transformed Gaussian

function x_out = Gaussian_mean_NMM(x,P,params)

N_params = params.N_params;
N_states = size(x,1) - N_params;
N_synapses = N_states/2;
N_inter_synapses = params.N_inter_synapses;
Tot_inter_synapses = sum(N_inter_synapses);

% pull all the parameters out of the struct
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A_mat = params.A_mat;
alpha = params.alpha;
tau = params.tau;
input_index = params.input_index;
connection_mat = params.connection_mat;
varsigma = params.varsigma;
v0 = params.v0;
external_input = params.input_mu;
external_in_scale = params.external_in_scale;
dt = params.dt;

% grab the membrane potentials and connectivity coefficients
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% alpha = ***
v_in = x(1:2:N_states-2*Tot_inter_synapses,:);
delayed_fr = x(N_states-2*(Tot_inter_synapses-1):2:N_states,:);

alpha = alpha';     % this is just to get the dimension right 
                    % for later when we augment
                    % the state vector

% calculate expected firing rates of each of the populations
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
intra_connection_mat = connection_mat(:,1:end-Tot_inter_synapses);
summed_membrane_voltage = zeros(size(intra_connection_mat,1),size(x,2));      % initialise for speed
for k=1:size(x,2)
    summed_membrane_voltage(:,k) = intra_connection_mat*v_in(:,k);
end


% need the variance of the mean membrane voltages
Pv = P(1:2*(N_states+1):end-2*Tot_inter_synapses);   % this is just the variance for v_in
voltage_var = connection_mat * Pv';
% this is the firing rate of each of the populations
[firing_rates,~] = exp_sigmoid_erf(varsigma,v0,summed_membrane_voltage,voltage_var);
% this is the index for the delayed firing rates (which region they are
% going to)
get_delayed_fr = cumsum([0 N_inter_synapses]);

% This is the linear part of the model
%
if N_params > 0
    % concatenate with I if there are parameters to estimate
    A_mat = blkdiag(A_mat,eye(N_params));
end

% This is the nonlinear part of the model
% 
gx = zeros(2*N_synapses+N_params,size(x,2));
% only put values in where the synapses are
region = 1;     % *** ASSUMES ONE EXTERNAL INPUT PER REGION. TRY AND FIX
for n = 1:N_synapses
    % the input is constant
    % NB the delayed firing rates are treated as additive constants with
    % the external input. The total input is scaled to be a proportion of
    % local and distant connections
    if input_index(n) == 0        
        if Tot_inter_synapses ~= 0
            gx(2*n,:) = 1./tau(n).*alpha(n,:).*dt.*((1-external_in_scale(region))*external_input(n) ...
                + external_in_scale(region)* ...
                sum(delayed_fr(get_delayed_fr(region)+1:get_delayed_fr(region+1))));
        else
            gx(2*n,:) = 1./tau(n).*alpha(n,:).*dt.*external_input(n);
        end
        region = region+1;
        
        % otherwise get the correct firing rate
        % the augmented states should be indexed in connection order
    else
        gx(2*n,:) = 1./tau(n).*alpha(n,:).*dt.*firing_rates(input_index(n),:);
    end
end

% calculate x_{n+1}
%
x_out = A_mat*x + gx;

end