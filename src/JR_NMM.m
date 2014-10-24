function [params, syn_index] = JR_NMM(params, new_params, syn_index, col)

P_max = params.P_max;
S_max = params.S_max;

alpha = params.alpha;
tau = params.tau;
constant_input = params.constant_input;
input_index = params.input_index;
connection_mat = params.connection_mat;
H = params.H;

%These values can be altered
He = new_params.He;
Hi = new_params.Hi;
taue = new_params.taue;
taui = new_params.taui;
taud = new_params.taud;
K = new_params.K;
noise_mean = new_params.noise_mean;

%These connection strengths aren't changed
JR_C1 = params.JR_C1;
JR_C2 = params.JR_C2;
JR_C3 = params.JR_C3;
JR_C4 = params.JR_C4;

e0 = params.e0;

% label the populations from 1 to N_populations
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% pop 2 = excitatory interneuron
% pop 1 = superficial pyramidal
% pop 3 = superficial inhibitory interneuron
% pop__ = external



%% synapse from external input to Superficial Py
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
syn_index = syn_index + 1;

constant_input(col+1) = noise_mean;

C(syn_index) = 1;                 % C1
gain = He;
tau(syn_index) = taue;
alpha(syn_index) = C(syn_index)*gain;           % note, there is no 2*e0 here

input_index(syn_index) = 0;             % from external input
population_to = col*P_max + 1;                          % to sPC
connection_mat(population_to,syn_index) = 1;            % now just make the connection

H(syn_index) = 1;                   % set to one if it contributes to the EEG (synapse to Py cells)

%% synapse from neighbouring column 1 to Superficial Py
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
syn_index = syn_index + 1;

C(syn_index) = K/2;
gain = He;
tau(syn_index) = taud;
alpha(syn_index) = C(syn_index)*gain;  

input_index(syn_index) = -1;             % from external input
population_to = col*P_max + 1;                          % to sPC
connection_mat(population_to,syn_index) = 1;            % now just make the connection

H(syn_index) = 1;                   % set to one if it contributes to the EEG (synapse to Py cells)

%% synapse from neighbouring column 2 to Superficial Py
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
syn_index = syn_index + 1;

C(syn_index) = K/2;
gain = He;
tau(syn_index) = taud;
alpha(syn_index) = C(syn_index)*gain;

input_index(syn_index) = -2;             % from external input
population_to = col*P_max + 1;                          % to sPC
connection_mat(population_to,syn_index) = 1;            % now just make the connection

H(syn_index) = 1;                   % set to one if it contributes to the EEG (synapse to Py cells)

%% synapse from excitatory interneurons (EIN) to superficial pyramidal cells (sPC)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
syn_index = syn_index + 1;

C(syn_index) = JR_C2;             % C2
gain = He;
tau(syn_index) = taue;
alpha(syn_index) = C(syn_index)*gain*2*e0;

input_index(syn_index) = col*P_max + 2;         % from EIN
population_to = col*P_max + 1;      % sPC
connection_mat(population_to,syn_index) = 1;

H(syn_index) = 1;                   % set to one if it contributes to the EEG (synapse to Py cells)

%% synapse from superficial pyramidal cells (sPC) to superficial inhibitory interneurons (sIIN)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
syn_index = syn_index + 1;

C(syn_index) = JR_C3;             % C3
gain = He;                  % excitatory gain
tau(syn_index) = taue;
alpha(syn_index) = C(syn_index)*gain*2*e0;

input_index(syn_index) = col*P_max + 1;         % from sPC
population_to = col*P_max + 3;      % to sIIN
connection_mat(population_to,syn_index) = 1;

H(syn_index) = 0;                   % set to one if it contributes to the EEG (synapse to Py cells)

%% synapse from superficial inhibitory interneurons (sIIN) to superficial pyramidal cells (sPC)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
syn_index = syn_index + 1;

C(syn_index) = JR_C4;             % C4
gain = Hi;              % inhibitory gain
tau(syn_index) = taui;
alpha(syn_index) = C(syn_index)*gain*2*e0;

input_index(syn_index) = col*P_max + 3;         % from sIIN
population_to = col*P_max + 1;      % sPC
connection_mat(population_to,syn_index) = 1;

H(syn_index) = 1;                   % set to one if it contributes to the EEG (synapse to Py cells)

%% superficial pyramidal cells (sPC) to excitatory interneuron (EIN)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
syn_index = syn_index + 1;

C(syn_index) = JR_C1;               % C15
gain = He;                  % excitatory gain
tau(syn_index) = taue;
alpha(syn_index) = C(syn_index)*gain*2*e0;

input_index(syn_index) = col*P_max + 1;         % from sPC
population_to = col*P_max + 2;      % to EII
connection_mat(population_to,syn_index) = 1;

H(syn_index) = 0;                   % set to one if it contributes to the EEG (synapse to Py cells)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%
% if syn_index ~= S_max
%     error('N_synapses set incorrectly')
% end

params.alpha = alpha;
params.tau = tau;
params.constant_input = constant_input;
params.input_index = input_index;
params.connection_mat = connection_mat;
params.H = H;

end