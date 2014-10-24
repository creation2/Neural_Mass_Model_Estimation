% model parameters
% ~~~~~~~~~~~~~~~~~

params.scale = 100;       % scales membrane voltages
params.N_params = 0;      % used in parameter estimation
params.mean_est = 0;      % used in calculation of analytic mean

% Model configuration
% ~~~~~~~~~~~~~~~~~~~~~~
N_cols = 4; %number of columns in the model
P_max = 3;  %maximum number of populations per column
S_max = 7;  %maximum number of synapses per column
params.N_cols = N_cols;
params.P_max = P_max;
params.S_max = S_max;

% external input
% ~~~~~~~~~~~~~~~~~
% NOISE
%
params.constant_input = zeros(1,N_cols);
%params.noise_mean = [220 90 90 90];
params.noise_mean = 220;
params.noise_var = 5.74;

% INTER-COLUMN
%
params.K = [47 39 27 43];
params.k = [0.5 0.2 0.1 0.4];
%        SEIZURE
% params.K = [1 100 100 100];
% params.k = [0 0.5 0.5 0.5];

% connectivity parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~
JR_C = 135;
params.JR_C1 = JR_C;
params.JR_C2 = 0.8*JR_C;
params.JR_C3 = 0.25*JR_C;
params.JR_C4 = 0.25*JR_C;

% sigmoid parameters
% ~~~~~~~~~~~~~~~~~~
params.v0 = 6e-3;              % threshold parameter
params.e0 = 2.5;               % half max firing rate - mid point of sigmoid     
params.varsigma = 0.003034; %560      % sigmoid slope erf/exp

% synapse parameters
% ~~~~~~~~~~~~~~~~~~
params.He = 3.25e-3;           % gain of excitatory synapses
params.Hi = -22e-3;            % gain of inhibitory synapses
params.taui = 0.020;           % inhibitory time constant
params.taue = 0.010;           % inhibitory time constant

params.taud = 0.0303;         % intercolumn delay time constant

% initalise column parameters
% ~~~~~~~~~~~~~~~~~
syn_index = 0;
params.alpha = zeros(1,N_cols*S_max);
params.tau = zeros(1,N_cols*S_max);
params.H = zeros(1,N_cols*S_max);
params.connection_mat = zeros(N_cols*P_max, N_cols*S_max);
params.input_index = zeros(1,N_cols*S_max);

