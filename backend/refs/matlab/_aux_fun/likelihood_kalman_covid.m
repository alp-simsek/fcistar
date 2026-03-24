function  loglik = likelihood_kalman_covid(y,x,theta,map_theta_to_mats,settings)

% Inputs

% y = T x n_y matrix of n_y measurment time series observations
% x = T x n_x matrix of n_x exogenous variables. This includes constants and lags.
% theta = n x 1 vector, includes all parameters.

%   State Equation
%   xi_{t+1} = F xi_{t} + v_{t+1}
%
%   Measurement Eq
%   y_t = A' x_t + H' xi_{t} + w_t
%
%   Error covariance
%   Q = E[v_t v_t']
%   R = E[w_t w_t']

%this function is application-specific, and maps the paramter vector theta to the state-space represtation
% given above.

fixed_param = settings.fixed_param;
covid_dummies = settings.covid_dummies;
matrices = map_theta_to_mats(theta,fixed_param,covid_dummies);

if settings.detrend_y == 1
    y(:,1) = y(:,1)-theta(5)*(1:length(y))';
    x(:,1) = x(:,1)-theta(5)*(0:length(y)-1)';
    x(:,2) = x(:,2)-theta(5)*(-1:length(y)-2)';
elseif settings.detrend_y == 2
    y(:,1) = y(:,1)-theta(5)*(1:length(y))';
    x(:,1) = x(:,1)-theta(5)*(0:length(y)-1)';
    x(:,2) = x(:,2)-theta(5)*(-1:length(y)-2)';
    x(:,3) = x(:,3)-theta(5)*(-2:length(y)-3)';
end


F = matrices.F;
A = matrices.A;
H = matrices.H;
Q = matrices.Q;
R = matrices.R;
kappa_seq = matrices.kappa_seq;

settings.smoother = 0;
% run the Kalman filter and obtain everything
output = kalman_filter_covid(y,x,F,A,H,Q,R,kappa_seq,settings);

% compute the log-likelihood.

loglik_vec = output.loglik_vec;
loglik = sum(loglik_vec); 
end

