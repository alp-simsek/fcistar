function V_mat_aux = compute_opg(theta_opt,constraints,evaluate_kalman)

% This script computes the standard error of parameters and states
% (accounting for parameter estimation error) using the OPG to get the
% variance-covariance matrix of the MLE estimator, and then using the
% monte carlo procedure in Hamilton (1986) and also used by HLW to compute
% the SEs for the states.

% INPUTS

% theta_opt N x 1 vector of the MLE estimate
% constraints N x 1 vector of logical values. If an element is = 1 it means
% the parameter hit the constraint and thus we don't report SEs for it.
% evaluate_kalman: function that computes the Kalman Filter. We are
% interested in the log-likelihood vector. This should be only a function
% of theta.

% First, compute the variance covariance matrix using OPG.
UB = constraints.UB;
LB = constraints.LB;
cons_vec = (theta_opt==LB)+(theta_opt==UB);

unconst_cols = find(cons_vec==0);
const_cols = find(cons_vec==1);
delta = 1e-5; %this is the step size for numerical derivatives.
out_base = evaluate_kalman(theta_opt);
lik_vec_base = out_base.loglik_vec;


T = length(lik_vec_base);
N_param = length(theta_opt);
lik_gradient = NaN(T,N_param);

for i = 1:N_param
if cons_vec(i)==0
    theta_prime = theta_opt;
    theta_prime(i) = theta_opt(i)+delta;
    output_prime = evaluate_kalman(theta_prime);
    lik_gradient(:,i) = (output_prime.loglik_vec-lik_vec_base)/delta;
end
end

% compute information matrix from the outer product of gradients

I_mat =  (lik_gradient'* lik_gradient)/T;
I_mat_aux = I_mat(unconst_cols,unconst_cols);
V_mat_aux = inv(I_mat(unconst_cols,unconst_cols))/T;
V_mat_aux_2 = (1/T)*(I_mat_aux\eye(length(unconst_cols)));
end



