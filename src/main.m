% column model

clc
close all
clearvars

% Turn randomness on/off
% ~~~~~~~~~~~~~~~~~~~~~~~~~
%
rng(0);

% Simulation Parameters
% ***********************

T = 5;                                  % seconds
Fs = 1e3;                               % sample rate (samples/second)
dt = 1/Fs;                              % this is the time step for integration
N_samples = T*Fs;                       % number of samples to generate
t = 0:dt:(N_samples-1)*dt;              % vector for all time points
% ONLY mV working for now
volts = 0;                              % 1 for V, 0 for mV


% this is where all the parameters of the neural masses are
initialise_params;                
% initalise for speed
x_out = zeros(Tot_states,N_samples);

% Simulation Noise
% *****************
% This is the external input noise
%
sigma = params.input_var;
input_index = params.input_index;
Q1 = zeros(size(x_out,1),1);
Q = zeros(Tot_synapses,1);
Q(input_index == 0) = (sqrt(dt)*sigma*He/taue).^2;
Q1(2:2:end) = Q;
Q1 = diag(Q1);

% This is the measurement nosie
%
R1 = 0.3;                   % measurement noise std in mV

% Generate trajectory
%
disturbance = mvnrnd(zeros(size(x_out,1),1),Q1,N_samples)';
for n=1:N_samples-1
    % integrate the model over time
    x_out(:,n+1) = neural_model(x_out(:,n),params) + disturbance(:,n);
end

% Make measurements
%
y = H * x_out;

if volts==1
    R = (1e-3*R1/(1/100))^2;
else
    %R = (R1)^2;
    R = 0;
end
obs_noise = sqrt(R)*randn(N_regions,N_samples);
%R = 0.009;
y = y + obs_noise;

figure;
if N_regions > 1
    for n = 1:N_regions
        subplot(2,2,n), plot(t,10*y(n,:),'k')
        box off
    end
else
    plot(t,10*y,'k')
    box off
end

%% Filter
N_states = Tot_states;
N_params = 0;
N_aug = N_states + N_params;

params2 = params;
params2.N_params = N_params;

% FILTER PARAMETERS
alpha = 1;                  % between 0 and 1. Use 1 for large N_aug
beta = 2;                   % 2 is optimal for Gaussians
kappa = 3 - N_states;       % ***Don't change***** 3 or 4
mat = 0;                    % zero

% initialise for speed
x_hat = zeros(N_states,N_samples);
P_hat = zeros(N_states,N_samples);

% inital guess of mean (M) and cov (P) for x
M = zeros(N_states,1);
P = diag(var(x_out,1,2));

x_hat(:,1) = M;
P_hat(:,1) = diag(P);

% this is for numerical stability
Q = Q1;
%Q(Q == 0) = 1e-14;

constraints = [repmat([0 0.1],4,1); ...   
    repmat([0 20],4,1); ...
    repmat([0 10],4,1); ...              
    repmat([-40 0],4,1); ...
    repmat([0 20],4,1); ...
    repmat([0 5],4,1); ...
    repmat([0 5],4,1)];

% Index of the states which are constrained (alpha parameters)

key = N_states+1:N_aug;

for n = 1:N_samples-1
    
   M = Gaussian_mean_NMM(M,P,params);
   [~,P] = ukf_predict1(M,P,@(x)neural_model(x,params2),Q,alpha,beta,kappa,constraints,key);
   [M,P] = KF_update(M,P,y(:,n+1),H,R);
   
   x_hat(:,n+1) = M;
   P_hat(:,n+1) = diag(P);
end

%% Plots
colours = ['m' 'b' 'r' 'k' 'g'];
% alpha = params.alpha;
% figure;
% for n = 1:4
%     plot(t,alpha(7*(n-1)+3)*ones(size(t)), colours(n), ...
%         t,x_hat(end-(4-n),:),colours(n));
%     hold on
% end
% figure;
% for n = 1:4
%     plot(t,alpha(7*(n-1)+2)*ones(size(t)), colours(n), ...
%         t,x_hat(end-(8-n),:),colours(n));
%     hold on
% end
for m = 1:4
    figure;
    for n = 1:5
        subplot(5,1,n),plot(t,x_out(2*n-1 + 10*(m-1),:),colours(n),t,x_hat(2*n-1 + 10*(m-1),:), [colours(n) '--']);
    end
end

figure;
for n = 1:8
    subplot(4,2,n),plot(t,x_out(end-n+1,:),colours(ceil(n/2)),t,x_hat(end-n+1,:), [colours(ceil(n/2)) '--']);
end