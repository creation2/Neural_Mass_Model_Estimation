
% this code propogates v_t and z_t of a synapse in the model though time.
% the integration is done using a simple Euler method.

% inputs
% ~~~~~~
% alpha is the synaptic gain
% tau is the time constant parameter
% v_t is the post-synaptic membrane potential at the previous time
% z_t is the derivative of v_t
% w_t is only used with habituation
% habituation = 0,1 to turn effect off/on
% mode = "transition" or "linear" for modelling or EKF estimation

% note: alpha, tau, v_t and z_t can be a scalar (for forward model) or
% vector (for sigma points in inverse model).


% outputs
% ~~~~~~~
% v is v_t after time dt
% z detivative of v
% w is updated habituation terms

function [z, v, w] = neural_mass(alpha,tau,input,v_t,z_t,w_t,habituation)

global dt

if habituation
    n1 = 20;                                            % depression rate
    n2 = 2;                                             % recovery rate
    v0 = 6e-3;                                          % half max firing rate - mid point of sigmoid
    varsigma = 560;                                     % sigmoid slope
    Q_max = 1 - 1./(1 + exp(varsigma*v0)); %4.8322;     % Qmax = the shifted maximum firing rate
    
    Anrr = 1 - w_t;
    if input > 0
        w = (-n1*input*w_t ./ Q_max + n2*Anrr)*dt + w_t;
    else
        w = (n2*Anrr)*dt + w_t;
    end
else
    w = w_t;
end

% the standard neural mass equation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
zeta = 1./tau;


v = (z_t*dt + v_t);
z = (w_t.*zeta.*alpha.*input - 2*zeta.*z_t - zeta.^2.*v_t)*dt + z_t;

