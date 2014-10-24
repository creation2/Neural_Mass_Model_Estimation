%UT_SIGMAS - Generate Constrained Sigma Points for Unscented Transformation
%
% Syntax:
%   X = ut_sigmas(M,P,c);
%
% In:
%   M - Initial state mean (Nx1 column vector)
%   P - Initial state covariance
%   c - Parameter returned by UT_WEIGHTS
%   constraints - matrix containing upper and lower bounds
%   key - logic vector indicating which states to apply constraints
%
% Out:
%   X - Matrix where 2N+1 sigma points are as columns
%
% Description:
%   Generates sigma points and associated weights for Gaussian
%   initial distribution N(M,P). For default values of parameters
%   alpha, beta and kappa see UT_WEIGHTS.
%
% See also UT_WEIGHTS UT_TRANSFORM UT_SIGMAS
%

% Copyright (C) 2006 Simo Särkkä
%
% $Id: ut_sigmas.m 109 2007-09-04 08:32:58Z jmjharti $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function X = ut_sigmas_constrained(M,P,c,constraints,key)

P = (P+P')/2;
A = chol(P + 1e-13*eye(length(M)))';
X = [zeros(size(M)) A -A];
X = sqrt(c)*X + repmat(M,1,size(X,2));

for n = 1:length(key)
    X_bound = X(key(n),:);
    X_bound(X_bound < constraints(n,1)) = constraints(n,1);
    X_bound(X_bound > constraints(n,2)) = constraints(n,2);
    X(key(n),:) = X_bound;
end

