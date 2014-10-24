clc
close all
clearvars -except rms_x rms_alpha bias_alpha window k

global dt

% Handle to model function
func_f  = @column_model_JR;

% Set up simulation variables
T = 100;                                % run time (seconds)
Fs = 1e3;                               % sample rate (samples/second)
dt = 1/Fs;%1e-3;                        % this is the time step for integration
N_samples = T*Fs;                       % number of samples to generate
t = 0:dt:(N_samples-1)*dt;              % vector for all time points

% Set up seizure transition
s_mode = 0;         % 1 for on, 0 for off
s_start = 30;       % start time (s)
s_length = 30;      % length of seizure (s)
s_transition = 5;   % time to transition (ramp params) into seizue (s)
s_gain = 2.5;       % increase He by this much

% Turn randomness on/off
% ~~~~~~~~~~~~~~~~~~~~~~~~
%
rng(0);

% Set up starting parameters
JRinitial_parameters;
% These params could be varied b/w regions
new_params.He = params.He;
new_params.Hi = params.Hi;
new_params.taue = params.taue;
new_params.taui = params.taui;
new_params.taud = params.taud;
new_params.noise_mean = params.noise_mean;

% Describe each column
for n = 1:N_cols
    % K is interconnectivity
    new_params.K = params.K(n);
    %new_params.noise_mean = params.noise_mean(n);
    % Now initialise each region
    % **************************
    [params, syn_index] = JR_NMM(params,new_params,syn_index,n-1);
end

N_states = 2*N_cols*S_max;              % total states
N_col_states = 2*S_max;                 % states per column
v = zeros(N_states/2,N_samples);        % initialise for speed - post-synaptic membrane potentials
z = zeros(N_states/2,N_samples);        % initialise for speed - derivatives of v
x_out = [v ; z];        %NB x_out actually alternating v1 z1 v2 z2 etc.


% Define the variance of the noise input
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Q = zeros(1,N_states/2);
sigma = params.noise_var;
He = params.He;             % gain of excitatory synapses
taue = params.taue;         % excitatory time constant

% Find the indices of the external input
% (index 0 is for synapses that receive u)
Q(params.input_index == 0) = (sqrt(dt)*sigma*He/taue)^2;
% The external input enters the membrane derivative equations
Q2 = zeros(1,N_states);
Q2(2:2:N_states) = Q;
Q = diag(Q2);
% additive random dissturbance
disturbance = mvnrnd(zeros(size(x_out,1),1),Q,N_samples)';

% Reshape the observation matrix
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N_obs = N_cols;                 % number of recording channels
H1 = zeros(1,N_states);
H1(1:2:N_states) = params.H;
H = zeros(N_obs,N_states);

if N_obs == N_cols
    for n = 1:N_cols
        ind = (n-1)*(N_col_states);
        ref = mod(n,N_cols) * N_col_states;
        H(n,ind+1:ind+N_col_states) = H1(ind+1:ind+N_col_states);
        H(n,ref+1:ref+N_col_states) = -H1(ref+1:ref+N_col_states);
    end
else
    H2 = H1(15:end); %***** magic no.
    for n = 1:N_cols-1
        ind = (n-1)*(N_col_states);
        ref = mod(n,N_obs) * N_col_states;
        H(n,ind+15:ind+2*N_col_states) = H2(ind+1:ind+N_col_states);
        H(n,ref+15:ref+2*N_col_states) = -H2(ref+1:ref+N_col_states);
    end
end

params.H = H;
params2 = params;

% Parameter Ramp for Seizure Generation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ***** PK THIS IS NOT GENERAL
if s_mode == 1
    % turn down the external input
    mu_ramp = 220:-26*dt:90+dt;
    %mu_ramp = repmat(mu_ramp,4,1);
    % turn up the synaptic gain (excitatory synapses only)
    alpha = params.alpha([1:5 7]);
    alpha_ramp = zeros(6,s_transition/dt);
    for n = 1:6
        alpha_ramp(n,:) = alpha(n):dt*(s_gain*alpha(n) - alpha(n))/s_transition:s_gain*alpha(n) - (alpha(n)/s_transition*dt);
    end
end

for n=1:N_samples-1                     % integrate the model over time
 
  
    % ******************* SEIZURE ********************************
    % PK THIS IS NOT GENERAL
    if s_mode == 1
        %RAMP INTO SEIZURE (col 1)
        if (n > s_start/dt) && (n <= (s_start+s_transition)/dt)
            params.alpha([1:5 7]) = alpha_ramp(:,n-s_start/dt)';
            %params.constant_input = mu_ramp(:,n-s_start/dt)';
            params.k(2:4) = 0.8;
            params.constant_input(1) = mu_ramp(n-s_start/dt);
        end
        %RAMP OUT
        s_end = (s_start+2*s_transition+s_length);
        if (n >= (s_start+s_transition+s_length)/dt) && (n < s_end/dt)
            params.alpha([1:5 7]) = alpha_ramp(:,s_end/dt-n)';
            %params.constant_input = mu_ramp(:,s_end/dt-n)';
            params.k(2:4) = 0.5; %0.4;
            params.constant_input(1) = mu_ramp(s_end/dt-n);
        end
    end
    % ************************************************************** 
    
    % Transition the model
    x_out(:,n+1) = func_f(x_out(:,n),params) + disturbance(:,n);
    
end

% make observation
y1 = H*x_out;

% add noise to measurements
R = .001;
obs_noise = sqrt(R).*randn(N_obs,N_samples); %sqrt(R).*randn(N_obs,N_samples);
y = y1 + obs_noise;

%% plot simulation
% figure;
% subplot(2,2,1), plot(t,10*y1(1,:));%,t,y(1,:),'rx');
% subplot(2,2,2), plot(t,10*y1(2,:));%,t,y(2,:),'rx');
% subplot(2,2,3), plot(t,10*y1(3,:));%,t,y(3,:),'rx');
% subplot(2,2,4), plot(t,10*y1(4,:));%,t,y(4,:),'rx');

%% Apply Kalman Filter
%~~~~~~~~~~~~~~~~~~~~~
%

N_params = 28;                  % no of parameters to estimate
params2.N_params = N_params;    % tell the filter how many parameters to estimate
N_aug = N_states + N_params;    % augmented state vector

% FILTER PARAMETERS
alpha = 1;                  % between 0 and 1. Use 1 for large N_aug
beta = 2;                   % 2 is optimal for Gaussians
kappa = 3 - N_states;       % ***Don't change***** 3 or 4
mat = 0;                    % zero

H = [H zeros(N_obs,N_params)];  % add parameters to observation function (unobserved)

x_hat = zeros(N_aug,N_samples);     % initialise augmented state vector
P_hat = zeros(N_aug,N_samples);     % initialise model covariance matrix (to store diagonals)

% first guess of P for the states is variance of the simulation
x_var = var(x_out,0,2);
col_var = x_var(1:2:end);
col_var = reshape(col_var,7,4)';
col_var = col_var ./ repmat(max(col_var),4,1);

P0 = [x_var; ...
    .1*ones(4,1).*col_var(:,1); ...
     10*ones(4,1).*col_var(:,4);...
     1*ones(4,1).*col_var(:,5);... 
     70*ones(4,1).*col_var(:,6);... %200
     60*ones(4,1).*col_var(:,7);... %150
     5*ones(4,1).*col_var(:,2);... 
     5*ones(4,1).*col_var(:,3)]; 

% P0 = [x_var; ...
%     .1*ones(4,1).*col_var(:,1); ...
%      80*ones(4,1).*col_var(:,4);...
%      1*ones(4,1);... 
%      200*ones(4,1);...
%      100*ones(4,1);...
%      .05*ones(4,1).*col_var(:,2);... 
%      .05*ones(4,1).*col_var(:,3)]; 

% % % SEIZURE SIMULATION
% P0 = [x_var; ...
%     .1*ones(4,1); ...
%     10*ones(4,1); ...
%     5*ones(4,1); ...
%     60*ones(4,1); ...
%     10*ones(4,1); ...
%     [.1;5;5;5]; ...
%     [.1;5;5;5]];  %.1 10 5 100 10 1 1

M = x_hat(:,1);   % first guess of states/params
P = diag(P0);     % first estimate of covariance

% Process noise, Q
estQ = zeros(N_aug,1);
estQ(1:N_states) = diag(Q);        % use the known process noise variance

% PARAMETER NOISE - FOR TRACKING
% **********************************
% estQ(N_states+1:N_states+4) = 1e-7*abs(params.alpha(1:7:end)); %1e-7,1e-5,1e-5...
% estQ(N_states+1+4:N_states+8) = 1e-5*abs(params.alpha(4:7:end));
% estQ(N_states+1+8:N_states+12) = 1e-5*abs(params.alpha(5:7:end));
% estQ(N_states+1+12:N_states+16) = 1e-5*abs(params.alpha(6:7:end));
% estQ(N_states+1+16:N_states+20) = 1e-5*abs(params.alpha(7:7:end));
% estQ(N_states+1+20:N_states+24) = 1e-5*abs(params.alpha(2:7:end));
% estQ(N_states+1+24:N_states+28) = 1e-5*abs(params.alpha(3:7:end));


estQ = diag(estQ);

% Apply bounds to augmented state vector
% *************************************
constraints = [repmat([0 0.1],4,1); ...   
    repmat([0 20],4,1); ...
    repmat([0 20],4,1); ...   %10           
    repmat([-40 0],4,1); ...
    repmat([0 20],4,1); ...
    repmat([0 5],4,1); ...
    repmat([0 5],4,1)];

% Index of the states which are constrained (alpha parameters)
%
key = N_states+1:N_aug;
%key = 0;

for n = 1:N_samples-1
      
    % evaluate variance for analytic mean
    state_var = get_variance(P(1:N_states,1:N_states),N_col_states,N_states);
    params2.state_var = state_var;   
   
    % ******************* SEIZURE ********************************
    if s_mode == 1
        %RAMP INTO SEIZURE (col 1)
        if (n > s_start/dt) && (n <= (s_start+s_transition)/dt)
            params2.alpha([1:5 7]) = alpha_ramp(:,n-s_start/dt)';
            params2.constant_input(1) = mu_ramp(n-s_start/dt);
            params.k(2:4) = 0.8;
        end
        %RAMP OUT
        s_end = (s_start+2*s_transition+s_length);
        if (n >= (s_start+s_transition+s_length)/dt) && (n < s_end/dt)
            params2.alpha([1:5 7]) = alpha_ramp(:,s_end/dt-n)';
            params2.constant_input(1) = mu_ramp(s_end/dt-n);
            params.k(2:4) = 0.5;
        end
    end
    % **************************************************************
    
    % Now update the parameter estimates
    params2 = update_params(x_hat(N_states+1:end,n)',[],[],params2);
    
    %apply filter
    [M,P] = ukf_predict1(M,P,func_f,estQ,params2,alpha,beta,kappa,mat,constraints,key);
    [M,P,K] = ukf_update1(M,P,y(:,n+1),H,R,[],alpha,beta,kappa,mat);
    
    %store next estimate and covariance
    x_hat(:,n) = M;
    P_hat(:,n) = diag(P);
    
end


%% Plots
colours = ['m' 'b' 'r' 'k'];
alpha = params.alpha;
figure;
for n = 1:4
    plot(t,alpha(7*(n-1)+3)*ones(size(t)), colours(n), ...
        t,x_hat(end-(4-n),:),colours(n));
    hold on
end
figure;
for n = 1:4
    plot(t,alpha(7*(n-1)+2)*ones(size(t)), colours(n), ...
        t,x_hat(end-(8-n),:),colours(n));
    hold on
end
% % % for m = 1:N_cols
% % %
% % %     figure;
% % %     for n = 1:7
% % %         subplot(7,1,n),plot(t,x_out(2*n-1 + 14*(m-1),:),colours(n),t,x_hat(2*n-1 + 14*(m-1),:), [colours(n) '--']);
% % %     end
% % % end