% specifics the initial configuration of the cortical circuit
%
% This code sets up vectors of parameters for the model. This is what
% will be moved into a GUI eventually to interactively set the values.

% PK At the moment it is set up for identical regions. If the regions are
% different there should probably be a separate function to set the gains
% and time constants and leave this script for just putting it all together

% In this code the populations are ordered SPy (1st) EIN (2nd) IIN (3rd)

params.N_params = 0;    %number of parameters to estimate (initally 0 to forward model)
if volts
    params.scale = 100;
else
    params.scale = 1;
end

N_regions = 4;
N_populations = [3 3 3 3];      % not including the external input
N_synapses = [5 5 5 5];         % intra-region

% This is where you define the connectivity between the models
% want an M by M matrix (for M regions) and put a 1 where the connection is
% made (unidirectional) if the matrix is symmetric you have bidirectional
% connections.
% This matrix will be made in GUI eventually
% connections will be made with a synapse that is Py-Py
% *** SHOULD CHANGE TO PY-EX
% There is only one type of connection so far. Could add others using the
% same thing with a different matrix maybe? Long range etc
if N_regions > 1
    region_connections = [0 1 0 1; ...
        1 0 1 0; ...
        0 1 0 1; ...
        1 0 1 0];
end

% Inter-region scale
% *******************

% The same length as number of regions

% K is the scale for the incoming connections to the region. It is divided
% by the number of inputs for the region
% k is the relative contribution of external input that is local or distant
% (1-k)*distant + k*local for 0 < k < 1
K = [47 39 27 43];
%K = [20 20 20 20];
k = [0.25 0.25 0.25 0.25];

Tot_populations = sum(N_populations);
Tot_synapses = sum(N_synapses);
Tot_states = 2*Tot_synapses;

% external input        % at the moment input is to Py. SHOULD CHANGE
% ***************
const_input = 300*ones(1,N_regions);        % firing rate in APs per second
const_var = 5.74*ones(1,N_regions);

% Here are the constant parameters for the model

% sigmoid parameters
% ~~~~~~~~~~~~~~~~~~
v0 = 6*ones(Tot_populations,1);             % threshold parameter
e0 = 2.5;                             % half max firing rate - mid point of sigmoid
r = 0.56*ones(Tot_populations,1);           % sigmoid slope

if volts
    v0 = 1e-3.*v0;
    r = 1e3.*r;
end

varsigma = 1.699./r;            % correct for the erf sigmoid

% JR scale factors
% ****************
JR_C = 135;
JR_C1 = JR_C;
JR_C2 = 0.8*JR_C;
JR_C3 = 0.25*JR_C;
JR_C4 = 0.25*JR_C;

% JR alpha function
% *****************
He = 3.25;           % gain of excitatory synapses
Hi = -22;            % gain of inhibitory synapses
taui = 0.020;        % inhibitory time constant
taue = 0.010;        % inhibitory time constant
taud = 0.0303;       % interregional delay time constant (for coupled models)
if volts
    He = 1e-3.*He;
    Hi = 1e-3.*Hi;
end

% Population Connectivity
%
% (py = 1, ex = 2, in = 3)
% this has the index of the pyramdial populations for every region
% py is always the first population in the region
go_to_next_pop = 1+cumsum([0 N_populations(1:end-1)])';

% Intra-regional
% the template is for a single region. This code creates identical coupled regions
% population each synapse comes from
from_index = zeros(1,Tot_synapses);
template_from = [0 2 1 3 1];  % set to 0 for external input, otherwise the population index
% population each synapse goes to (pouplation number)
to_index = zeros(1,Tot_synapses);
template_to = [1 1 3 1 2];

% This loop here just repeats the connectivity template for every region
% will need to be changed when the regions are not identical
from_index(1:N_synapses(1)) = template_from;
to_index(1:N_synapses(1)) = template_to;
for n = 1:N_regions-1
    ind = sum(N_synapses(1:n));
    next_from_index = sum(N_populations(1:n))+template_from;
    next_from_index(next_from_index == sum(N_populations(1:n))) = 0;
    from_index(ind+1:ind+N_synapses(n+1)) = ...
        next_from_index;
    
    to_index(ind+1:ind + N_synapses(n+1)) = ...
        sum(N_populations(1:n))+template_to;
end

% Inter-regional
%
if N_regions > 1
    % here add the new synapses to update the model
    N_inter_syn = sum(region_connections);
    Tot_inter_syn = sum(N_inter_syn);
    N_synapses = N_synapses + N_inter_syn;
    Tot_synapses = Tot_synapses + Tot_inter_syn;
    Tot_states = 2*Tot_synapses;
    
    % find which regions they connect
    [region_to,region_from] = find(region_connections);
    % make sure the index is for Py instead of region number
    inter_to_index = go_to_next_pop(region_to)';
    inter_from_index = go_to_next_pop(region_from)';
    % now add the new synapses to the others
    % keep them on the end because they are actually delayed firing rates so
    % need to process them separately to the other states
    from_index = [from_index inter_from_index];
    to_index = [to_index inter_to_index];
else
    N_inter_syn = [];
    Tot_inter_syn = 0;
end

% Synapse parameters
%
% Here we make the synapses distinct using alpha and tau
% This bit also assumes identical regions
% time constants
% **************
tau = [taue taue taue taui taue];
tau = repmat(tau,1,N_regions);
tau = [tau taud*ones(1,Tot_inter_syn)];      % the inter-region synapses all on the end

% connectivity gains
% *************************
alpha_intra = [He ...  % u(t) to SPy
    2*e0*JR_C2*He ...  % EIN to SPY
    2*e0*JR_C3*He ...  % SPy to IIN
    2*e0*JR_C4*Hi ...  % IIN to SPy
    2*e0*JR_C1*He];    % SPy to EIN
alpha_intra = repmat(alpha_intra,1,N_regions);
alpha_inter = zeros(1,Tot_inter_syn);

if N_regions > 1
    ind = 0;
    for n = 1:N_regions
        % need to know how many inputs for the region
        num_inputs = sum(region_to == n);
        alpha_inter(ind+1:ind+num_inputs) = (He*K(n)/num_inputs) * ones(1,num_inputs);
        ind = ind+num_inputs;
    end
end

alpha = [alpha_intra alpha_inter];

% Observation Function
%
% we only observe pyramidal populations (which will be indexed as the first pop
% for each region)
H = zeros(N_regions,Tot_states);
H1 = repmat(to_index,N_regions,1) == repmat(go_to_next_pop,1,Tot_synapses);
H(:,1:2:end) = H1;        % NB the derivative states aren't measured
% NB we won't include the delay states in the measurement
H(:,end-2*Tot_inter_syn+1:end) = 0;

%% Initialise the model
% *********************
% external input
% only apply external input to the right synapses
external_input = zeros(1,Tot_synapses);
external_input(from_index == 0) = const_input;
% build up adjacency matrix
%
connection_mat = zeros(Tot_populations,Tot_synapses);
for n = 1:Tot_populations
    connection_mat(n,:) = to_index == n;
end

% Make the linear part of the model
% (this doesn't change so it is faster to initialize once)
A_mat = zeros(Tot_states);
B_mat = zeros(Tot_states);
zeta = 1./tau;
for n = 1:Tot_synapses
    A_mat(2*n-1:2*n,2*n-1:2*n) = [1 dt;
        -dt*zeta(n)^2 1-2*zeta(n)*dt];
    
    B_mat(2*n,2*n) = dt*alpha(n)*zeta(n);
end

params.A_mat = A_mat;
params.B_mat = B_mat;
params.alpha = alpha;
params.tau = tau;
params.input_index = from_index;
params.connection_mat = connection_mat;
params.N_inter_synapses = N_inter_syn;
params.external_in_scale = k;
params.varsigma = varsigma;
params.v0 = v0;
params.input_mu = external_input;
params.input_var = const_var;
params.dt = dt;