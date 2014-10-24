%Models a neural mass column/s
%
% ************ INPUTS *************
% x - state vector of membrane voltages and derivatives (alternating
% i.e. [v1 z1 ... vN zN]
% x_delay - state vector delayed by prop_delay to use as inter-column
% inputs
% params - struct containing model parameters

%************* OUTPUTS ***************
% x_out - next state vector

function out = column_model_JR(x, params)

N_params = params.N_params;
expectation = params.mean_est;

N_states = size(x,1) - N_params;
N_synapses = N_states/2;

%if there are parameters to be estimated
if N_params
    %get the parameter estimates
    est_params = x(N_states+1:end)';
    %use the estimates to update the struct
    params = update_params(est_params,N_params,N_synapses,params);
    %get the votlages & derivatives
    x = x(1:N_states);
end

% pull all the paramters out of the struct
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
P_max = params.P_max;
S_max = params.S_max;
N_cols = params.N_cols;

% connectivity parameters
connection_mat = params.connection_mat;
input_index = params.input_index;
k = params.k;
scale = params.scale;

% column parameters
alpha = params.alpha;
tau = params.tau;
varsigma = params.varsigma;
v0 = params.v0;


% grab the membrane potentials and the derivatives out of the struct
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% if mode == 1                            % forward model with additive noise on states / no habituation
%     v_in = x(1:N_synapses,:);
%     z_in = x(N_synapses+1:end,:);
%     w_in = ones(1,N_synapses);
%     habituation = zeros(1,N_synapses);  % turn off the habituation effect
% end
%
% if mode == 2                            % inverse model with additive noise on states / no habituation
%     v_in = x(1:N_synapses,:);
%     z_in = x(N_synapses+1:2*N_synapses,:);
%     w_in = x(2*N_synapses+1:end,:);
%     habituation = zeros(1,N_synapses);  % turn off the habituation effect
% end

%find the voltages in the state vector
v_in = x(1:2:N_states)./scale;
%Then the derivatives
z_in = x(2:2:N_states);
%habituation
w_in = ones(1,N_synapses);
habituation = zeros(1,N_synapses);  % turn off the habituation effect

%This is the mean membrane pot of each population
summed_membrane_voltage = connection_mat*v_in;


% run each synapse seperately
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
z_out = zeros(N_synapses,1);                                    % initialise for speed
v_out = zeros(N_synapses,1);
w_out = zeros(N_synapses,1);

%Get the firing rates from the membrane potentials
firing_rates = sigmoid_JR(varsigma,v0,summed_membrane_voltage);

if expectation
    %**************** CALCULATIONS FOR ANALYTIC MEAN ******************************
    %get the variance
    state_var = params.state_var;
    
    for n = 1:N_synapses
        %work out which column is contributing input
        current_col = ceil(n/S_max);
        
        if input_index(n) == 0
            input = (1-k(current_col))*params.constant_input(current_col);
            
        elseif input_index(n) == -1
            next_pop = mod(P_max*current_col + 1 , P_max*N_cols);
            input = mean_sigmoid_JR(varsigma,v0,summed_membrane_voltage(next_pop), ...
                state_var(next_pop));
            
        elseif input_index(n) == -2
            prev_pop = mod(P_max*current_col-5, P_max*N_cols);
            input = mean_sigmoid_JR(varsigma,v0,summed_membrane_voltage(prev_pop), ...
                state_var(prev_pop));
        else
            ind = input_index(n);
            input = mean_sigmoid_JR(varsigma,v0, ...
                summed_membrane_voltage(ind),state_var(ind));
        end
        [z_out(n), v_out(n), w_out(n)] = neural_mass(alpha(n),tau(n), ...
            input,v_in(n),z_in(n),w_in(n),habituation(n));
    end
    
else
    %**************** CALCULATIONS FOR FORWARD MODEL ************************
    
    for n = 1:N_synapses
        %work out which column is contributing input
        current_col = ceil(n/S_max);
        %external noise input
        if input_index(n) == 0
            input = (1-k(current_col))*params.constant_input(current_col);
            %input from adjacent col
        elseif input_index(n) == -1
            next_pop = mod(P_max*current_col + 1 , P_max*N_cols);
            input = firing_rates(next_pop);
            %input form adjacent col
        elseif input_index(n) == -2
            prev_pop = mod(P_max*current_col-5, P_max*N_cols);
            input = firing_rates(prev_pop);
        else
            input = firing_rates(input_index(n));
        end
        [z_out(n), v_out(n), w_out(n)] = neural_mass(alpha(n),tau(n), ...
            input,v_in(n),z_in(n),w_in(n),habituation(n));
    end
    
end

%     if mode == 1
%         z_out = z_out + sqrt(var_z(n))*randn;%noise(n);           % add noise to the appropriate states to simulate modelling error
%     end

% reconstruct state vector
% ~~~~~~~~~~~~~~~~~~~~~~~~

v_out = scale.*v_out;
out = zeros(N_states,1);
out(1:2:N_states) = v_out;
out(2:2:N_states) = z_out;
if N_params
    out = [out;est_params'];
end

end